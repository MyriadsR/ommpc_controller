#pragma once

#include "polynomial_trajectory.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

// Reconstruct a C3 trajectory from sampled position/velocity. Local polynomial
// fitting suppresses text rounding noise before differentiating to acceleration
// and jerk; Hermite segments preserve those derivatives between samples.
class TextTrajectory
{
public:
  struct Sample
  {
    Eigen::Vector3d p, v, a, j;
    double yaw;
  };

  bool load(const std::string &filename, double dt)
  {
    trajectory_.clear();
    yaws_.clear();
    if (!std::isfinite(dt) || dt <= 0.0) return false;
    std::ifstream file(filename);
    if (!file) return false;
    std::vector<Eigen::Vector3d> positions, velocities, accelerations, jerks;
    bool explicit_derivatives = false;
    std::string line;
    while (std::getline(file, line))
    {
      if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
      std::istringstream stream(line);
      std::vector<double> row;
      double value;
      while (stream >> value)
      { if (!std::isfinite(value)) return false; row.push_back(value); }
      if (!stream.eof() || (row.size() != 8 && row.size() != 14)) return false;
      if (positions.empty()) explicit_derivatives = row.size() == 14;
      if (explicit_derivatives != (row.size() == 14)) return false;
      if (explicit_derivatives)
      {
        accelerations.emplace_back(row[8], row[9], row[10]);
        jerks.emplace_back(row[11], row[12], row[13]);
      }
      positions.emplace_back(row[0], row[1], row[2]);
      velocities.emplace_back(row[3], row[4], row[5]);
      yaws_.push_back(row[6]);
    }
    if (positions.size() < 4) return false;
    dt_ = dt;
    const int count = positions.size();
    const int window = std::min(21, count);
    const int degree = std::min(7, 2 * window - 1);
    const double scale = (window - 1) * dt_ / 2.0;
    std::vector<Sample> knots(count);
    // The design depends only on the sample's location inside the window.
    std::vector<Eigen::MatrixXd> maps(window);
    for (int center = 0; center < window; ++center)
    {
      Eigen::MatrixXd design = Eigen::MatrixXd::Zero(2 * window, degree + 1);
      for (int k = 0; k < window; ++k)
      {
        const double x = (k - center) * dt_ / scale;
        for (int power = 0; power <= degree; ++power)
        {
          design(k, power) = std::pow(x, power);
          if (power > 0)
            design(window + k, power) = power * std::pow(x, power - 1);
        }
      }
      maps[center] = design.colPivHouseholderQr().solve(
          Eigen::MatrixXd::Identity(2 * window, 2 * window));
    }
    for (int i = 0; i < count; ++i)
    {
      if (explicit_derivatives)
      {
        // A 14-column flatness reference already supplies consistent a/j;
        // preserve them instead of differentiating rounded samples again.
        knots[i].p = positions[i]; knots[i].v = velocities[i];
        knots[i].a = accelerations[i]; knots[i].j = jerks[i];
        continue;
      }
      const int first = std::max(0, std::min(i - window / 2, count - window));
      Eigen::MatrixXd values(2 * window, 3);
      for (int k = 0; k < window; ++k)
      {
        values.row(k) = positions[first + k].transpose();
        values.row(window + k) = (velocities[first + k] * scale).transpose();
      }
      const Eigen::MatrixXd coeff = maps[i - first] * values;
      knots[i].p = coeff.row(0).transpose();
      knots[i].v = coeff.row(1).transpose() / scale;
      knots[i].a = 2.0 * coeff.row(2).transpose() / (scale * scale);
      knots[i].j = 6.0 * coeff.row(3).transpose() / (scale * scale * scale);
    }
    Eigen::Matrix4d end_matrix;
    for (int k = 0; k < 4; ++k)
    {
      const double power = k + 4;
      end_matrix.col(k) << 1.0, power, power * (power - 1),
          power * (power - 1) * (power - 2);
    }
    const Eigen::Matrix4d end_inverse = end_matrix.inverse();
    const int stride = explicit_derivatives ? 1 : std::min(10, count - 1);
    piece_dt_ = stride * dt_;
    for (int i = 0; i + 1 < count; i += stride)
    {
      const int next = std::min(i + stride, count - 1);
      const double segment_dt = (next - i) * dt_;
      Eigen::Matrix<double, 8, 3> coeff;
      coeff.row(0) = knots[i].p.transpose();
      coeff.row(1) = (knots[i].v * segment_dt).transpose();
      coeff.row(2) = (knots[i].a * segment_dt * segment_dt / 2.0).transpose();
      coeff.row(3) = (knots[i].j * segment_dt * segment_dt * segment_dt / 6.0).transpose();
      Eigen::Matrix<double, 4, 3> remainder;
      remainder.row(0) = knots[next].p.transpose() - coeff.topRows<4>().colwise().sum();
      remainder.row(1) = (knots[next].v * segment_dt).transpose()
          - coeff.row(1) - 2 * coeff.row(2) - 3 * coeff.row(3);
      remainder.row(2) = (knots[next].a * segment_dt * segment_dt).transpose()
          - 2 * coeff.row(2) - 6 * coeff.row(3);
      remainder.row(3) = (knots[next].j * segment_dt * segment_dt * segment_dt).transpose() - 6 * coeff.row(3);
      coeff.bottomRows<4>() = end_inverse * remainder;
      Eigen::MatrixXd piece(3, 8);
      for (int k = 0; k < 8; ++k)
        piece.col(7 - k) = coeff.row(k).transpose() / std::pow(segment_dt, k);
      trajectory_.emplace_back(segment_dt, piece);
    }
    return true;
  }

  double duration() const { return empty() ? 0.0 : dt_ * (yaws_.size() - 1); }
  bool empty() const { return trajectory_.getPieceNum() == 0; }

  Sample sample(double t) const
  {
    if (empty() || !std::isfinite(t))
      throw std::invalid_argument("Cannot sample an empty trajectory or non-finite time");
    Sample s;
    const double clamped = std::max(0.0, std::min(t, duration()));
    const int first = std::min(static_cast<int>(clamped / dt_), static_cast<int>(yaws_.size()) - 2);
    const int piece_index = std::min(static_cast<int>(clamped / piece_dt_), trajectory_.getPieceNum() - 1);
    const Piece &piece = trajectory_[piece_index];
    const double local = clamped - piece_index * piece_dt_;
    s.p = piece.getPos(local); s.v = piece.getVel(local);
    s.a = piece.getAcc(local); s.j = piece.getJer(local);
    const double alpha = std::max(0.0, std::min(1.0, clamped / dt_ - first));
    const double difference = std::atan2(std::sin(yaws_[first + 1] - yaws_[first]),
                                        std::cos(yaws_[first + 1] - yaws_[first]));
    s.yaw = yaws_[first] + alpha * difference;
    // Hold the final position after playback, without carrying the last
    // nonzero velocity/acceleration into an indefinitely moving reference.
    if (t >= duration()) { s.v.setZero(); s.a.setZero(); s.j.setZero(); }
    return s;
  }

private:
  Trajectory trajectory_;
  std::vector<double> yaws_;
  double dt_ = 0.01, piece_dt_ = 0.1;
};
