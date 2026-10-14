#include "ommpc_controller.hpp"
#include <gtest/gtest.h>
Parameter_t testParams()
{
  Parameter_t p{};p.step_T=.01;p.hover_percent=.7;
  p.Q_pos_x=p.Q_pos_y=p.Q_pos_z=1000;
  p.Q_velocity=20;p.Q_attitude_roll=p.Q_attitude_pitch=p.Q_attitude_yaw=40;
  p.R_thrust=p.R_roll=p.R_pitch=p.R_yaw=1;
  p.min_thrust=1;p.max_thrust=39;
  p.max_bodyrate_x=p.max_bodyrate_y=14;p.max_bodyrate_z=4;return p;
}
TEST(ThrustModel, LegacyLinearMapping)
{
  MpcController c;auto p=testParams();c.init(p);
  EXPECT_NEAR(c.thrustToAcceleration(.7),9.81,1e-12);
  EXPECT_NEAR(c.accelerationToThrust(9.81),.7,1e-12);
}
TEST(ThrustModel, QuadraticInverseAndAppliedBounds)
{
  MpcController c;auto p=testParams();
  p.thrust_model_power=2;p.thrust_model_offset=.1;p.thrust_model_gain=36.17032258064516;
  p.min_normalized_thrust=.02;p.max_normalized_thrust=.95;c.init(p);
  for(double u=.02;u<.96;u+=.01)
    EXPECT_NEAR(c.accelerationToThrust(c.thrustToAcceleration(u)),u,1e-12);
  c.setHoverReference(Eigen::Vector4d(0,0,4,0));
  EXPECT_NEAR(9.81-c.u_ub[0](0),1,1e-12);
  EXPECT_NEAR(9.81-c.u_lb[0](0),std::min(p.max_thrust,c.thrustToAcceleration(.95)),1e-12);
}
TEST(ThrustModel, ManeuverReserveDoesNotAffectHoverOrLanding)
{
  MpcController c;auto p=testParams();p.min_tracking_thrust=5;c.init(p);
  c.setHoverReference(Eigen::Vector4d(0,0,4,0));
  const double hover_min=9.81-c.u_ub[0](0);
  Odom_Data_t odom;
  std::vector<Eigen::Vector3d> positions(nstep+1,Eigen::Vector3d(0,0,4));
  std::vector<Eigen::Vector3d> zeros(nstep+1,Eigen::Vector3d::Zero());
  std::vector<double> yaw(nstep+1,0);
  c.setTextReference(positions,zeros,zeros,zeros,odom,0,yaw,true);
  EXPECT_NEAR(9.81-c.u_ub[0](0),5,1e-12);
  c.setTextReference(positions,zeros,zeros,zeros,odom,0,yaw);
  EXPECT_NEAR(9.81-c.u_ub[0](0),hover_min,1e-12);
}
int main(int argc,char**argv)
{
  ros::init(argc,argv,"thrust_model_test",ros::init_options::AnonymousName);
  ros::Time::init();
  ::testing::InitGoogleTest(&argc,argv);return RUN_ALL_TESTS();
}
