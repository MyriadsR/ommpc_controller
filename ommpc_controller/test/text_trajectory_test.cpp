#include "text_trajectory.h"
#include <gtest/gtest.h>
#include <cstdio>
#include <iomanip>
#include <unistd.h>

class TextTrajectoryTest : public ::testing::Test
{
protected:
  std::string filename;
  void SetUp() override
  {
    char name[] = "/tmp/ommpc_text_test_XXXXXX";
    const int descriptor = mkstemp(name);
    ASSERT_GE(descriptor, 0);
    close(descriptor);
    filename = name;
  }
  void TearDown() override { std::remove(filename.c_str()); }
  void writeCircle()
  {
    std::ofstream file(filename);
    file << std::fixed << std::setprecision(6);
    for (int i = 0; i <= 250; ++i)
    {
      const double t = i * .01, v = .5 + .125 * t;
      const double theta = (.5 * t + .0625 * t * t) / 1.2;
      file << 1.2 - 1.2 * cos(theta) << " 0 " << 3 + 1.2 * sin(theta)
           << " " << v * sin(theta) << " 0 " << v * cos(theta) << " 0 0\n";
    }
    file << "\n";
  }
};

TEST_F(TextTrajectoryTest, RoundedSamplesRecoverDynamicsAndConsistentDerivatives)
{
  writeCircle();
  TextTrajectory trajectory;
  ASSERT_TRUE(trajectory.load(filename, .01));
  EXPECT_NEAR(trajectory.duration(), 2.5, 1e-12);
  for (double t = .12; t < 2.38; t += .0173)
  {
    auto s = trajectory.sample(t);
    const double v = .5 + .125 * t, theta = (.5 * t + .0625 * t * t) / 1.2;
    Eigen::Vector3d acc(v * v / 1.2 * cos(theta) + .125 * sin(theta), 0,
                       -v * v / 1.2 * sin(theta) + .125 * cos(theta));
    Eigen::Vector3d jerk(3 * v * .125 / 1.2 * cos(theta) - pow(v, 3) / 1.44 * sin(theta), 0,
                        -3 * v * .125 / 1.2 * sin(theta) - pow(v, 3) / 1.44 * cos(theta));
    EXPECT_LT((s.a - acc).norm(), .002);
    EXPECT_LT((s.j - jerk).norm(), .05);
    const double epsilon = 1e-5;
    const auto left = trajectory.sample(t - epsilon), right = trajectory.sample(t + epsilon);
    EXPECT_LT(((right.p - left.p) / (2 * epsilon) - s.v).norm(), 1e-4);
    EXPECT_LT(((right.v - left.v) / (2 * epsilon) - s.a).norm(), 1e-4);
    EXPECT_LT(((right.a - left.a) / (2 * epsilon) - s.j).norm(), 1e-4);
  }
}

TEST_F(TextTrajectoryTest, DerivativesAreContinuousAcrossSegmentsAndEndpointHolds)
{
  writeCircle();
  TextTrajectory trajectory;
  ASSERT_TRUE(trajectory.load(filename, .01));
  for (double t = .1; t < 2.4; t += .1)
  {
    auto left = trajectory.sample(t - 1e-8), right = trajectory.sample(t + 1e-8);
    EXPECT_LT((left.p - right.p).norm(), 1e-6);
    EXPECT_LT((left.v - right.v).norm(), 1e-6);
    EXPECT_LT((left.a - right.a).norm(), 1e-6);
    EXPECT_LT((left.j - right.j).norm(), 1e-5);
  }
  auto end = trajectory.sample(3.0);
  EXPECT_TRUE(end.v.isZero()); EXPECT_TRUE(end.a.isZero()); EXPECT_TRUE(end.j.isZero());
}

TEST_F(TextTrajectoryTest, InvalidRowsAndSamplingPeriodAreRejected)
{
  std::ofstream file(filename); file << "1 2\n"; file.close();
  TextTrajectory trajectory;
  EXPECT_FALSE(trajectory.load(filename, .01));
  EXPECT_TRUE(trajectory.empty());
  writeCircle();
  EXPECT_FALSE(trajectory.load(filename, 0));
  EXPECT_FALSE(trajectory.load(filename, -.01));
}

TEST_F(TextTrajectoryTest, ExplicitDerivativesPreserveAnExactPolynomial)
{
  std::ofstream file(filename);
  file << std::fixed << std::setprecision(12);
  for(int i=0;i<=100;++i)
  {
    const double t=.01*i;
    file << t*t*t << " 0 3 " << 3*t*t << " 0 0 0 0 " << 6*t << " 0 0 6 0 0\n";
  }
  file.close();
  TextTrajectory trajectory;
  ASSERT_TRUE(trajectory.load(filename,.01));
  for(double t=.001;t<.99;t+=.0031)
  {
    auto s=trajectory.sample(t);
    EXPECT_NEAR(s.p.x(),t*t*t,1e-9);
    EXPECT_NEAR(s.v.x(),3*t*t,1e-8);
    EXPECT_NEAR(s.a.x(),6*t,1e-7);
    EXPECT_NEAR(s.j.x(),6,1e-5);
  }
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
