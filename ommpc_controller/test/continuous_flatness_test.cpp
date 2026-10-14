#include "continuous_flatness.h"
#include "so3_math.h"
#include <gtest/gtest.h>

TEST(ContinuousFlatness, FullPitchRotationsRemainAlignedAndContinuous)
{
  Eigen::Quaterniond previous = Eigen::Quaterniond::Identity();
  for (int i=0; i<=7200; ++i)
  {
    const double theta=i*M_PI/1800;
    const Eigen::Vector3d force(10*sin(theta),0,10*cos(theta));
    const Eigen::Vector3d jerk(40*cos(theta),0,-40*sin(theta));
    Eigen::Quaterniond q; Eigen::Vector3d omega;
    ASSERT_TRUE(continuousFlatnessReference(force,jerk,previous,q,omega));
    EXPECT_LT((q*Eigen::Vector3d::UnitZ()-force.normalized()).norm(),1e-8);
    EXPECT_GT(q.dot(previous),.999);
    EXPECT_LT((omega-Eigen::Vector3d(0,4,0)).norm(),1e-8);
    EXPECT_NEAR(q.norm(),1,1e-12);
    previous=q;
  }
}

TEST(ContinuousFlatness, AngularRateMatchesTransportOnThreeDimensionalPath)
{
  Eigen::Quaterniond previous=Eigen::Quaterniond::Identity();
  Eigen::Vector3d last_rate=Eigen::Vector3d::Zero();
  const double dt=1e-4;
  for(int i=0;i<50000;++i)
  {
    const double t=i*dt;
    Eigen::Vector3d force(sin(t),.2*sin(2*t),cos(t));
    Eigen::Vector3d jerk(cos(t),.4*cos(2*t),-sin(t));
    Eigen::Quaterniond q;Eigen::Vector3d rate;
    ASSERT_TRUE(continuousFlatnessReference(force,jerk,previous,q,rate));
    if(i>0)
    {
      const Eigen::Vector3d actual=SO3::log(previous.toRotationMatrix().transpose()*q.toRotationMatrix())/dt;
      EXPECT_LT((actual-.5*(last_rate+rate)).norm(),.001);
    }
    previous=q;last_rate=rate;
  }
}

TEST(ContinuousFlatness, AntipodalDirectionAndZeroThrustAreHandled)
{
  Eigen::Quaterniond q;Eigen::Vector3d w;
  ASSERT_TRUE(continuousFlatnessReference(-Eigen::Vector3d::UnitZ(),Eigen::Vector3d::Zero(),Eigen::Quaterniond::Identity(),q,w));
  EXPECT_LT((q*Eigen::Vector3d::UnitZ()+Eigen::Vector3d::UnitZ()).norm(),1e-10);
  EXPECT_FALSE(continuousFlatnessReference(Eigen::Vector3d::Zero(),Eigen::Vector3d::Zero(),q,q,w));
}

int main(int argc,char**argv) { ::testing::InitGoogleTest(&argc,argv);return RUN_ALL_TESTS(); }
