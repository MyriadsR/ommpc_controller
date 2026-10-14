#pragma once
#include <Eigen/Eigen>
#include <cmath>

// Lift a path of thrust directions continuously to SO(3). Transport the
// transverse axes instead of choosing a singular global yaw chart.
inline bool continuousFlatnessReference(const Eigen::Vector3d &force,
                                        const Eigen::Vector3d &jerk,
                                        const Eigen::Quaterniond &previous,
                                        Eigen::Quaterniond &attitude,
                                        Eigen::Vector3d &body_rate)
{
  // Zero collective force has no defined thrust axis; propagate failure to the
  // caller rather than inventing an attitude or dividing by a near-zero norm.
  if (!force.allFinite() || !jerk.allFinite() || force.norm() < 1e-6) return false;
  const Eigen::Vector3d z = force.normalized();
  const Eigen::Quaterniond seed = previous.normalized();
  const Eigen::Vector3d old_z = seed * Eigen::Vector3d::UnitZ();
  Eigen::Quaterniond transport;
  // At an antipodal jump, the shortest rotation axis is ambiguous. The old
  // reference body X axis supplies a deterministic transverse axis.
  if (old_z.dot(z) < -1.0 + 1e-10)
    transport = Eigen::Quaterniond(Eigen::AngleAxisd(M_PI, seed * Eigen::Vector3d::UnitX()));
  else
    transport = Eigen::Quaterniond::FromTwoVectors(old_z, z);
  attitude = (transport * seed).normalized();
  // q and -q are the same SO(3) attitude; keep one quaternion hemisphere for
  // adjacent reference samples so interpolation does not introduce a jump.
  if (attitude.dot(seed) < 0.0) attitude.coeffs() *= -1.0;
  const Eigen::Vector3d z_dot = (jerk - z * z.dot(jerk)) / force.norm();
  // Zero spin around z: omega_world = z x dz/dt. Convert that feedforward to
  // body coordinates, matching PX4's angular-rate command convention.
  body_rate = attitude.conjugate() * z.cross(z_dot);
  return attitude.coeffs().allFinite() && body_rate.allFinite();
}
