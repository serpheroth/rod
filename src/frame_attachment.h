#ifndef FRAME_ATTACHMENT_H
#define FRAME_ATTACHMENT_H
#include <functional>
#include "global.h"
#include "pbd_elastic_rod.h"
class FrameAttachment
{
public:
  typedef PBDElasticRod::Vec3 Vec3;
  typedef PBDElasticRod::Mat3 Mat3;
  typedef std::function<void (FrameAttachment&)> FrameUpdator;
  /**
   * @brief FrameAttachment frame components are the row vectors of the matrix
   * @param initial_frame
   * @param attached_frame_ frame of attached rod
   * @param pos postion of attached point
   */
  FrameAttachment(PBDElasticRod* rod, Mat3& initial_frame, Mat3& attached_frame, Vec3& pos);
  int UpdateFrame(FrameUpdator updator);
  void EnforceFrameAttachment(Vec3 &p0, Vec3 &p1, Vec3 &pe, real length);
  void EnforcePosAttachment(Vec3& p);

  inline void ComputeFrameDerivative(Vec3 &p0, Vec3 &p1, Vec3 &p2, Mat3& d,
                                     Mat3& d1p0, Mat3& d1p1, Mat3& d1p2,
                                     Mat3& d2p0, Mat3& d2p1, Mat3& d2p2,
                                     Mat3& d3p0, Mat3& d3p1, Mat3& d3p2);
  PBDElasticRod* rod_;
  Mat3 initial_frame_;
  Mat3 rod_frame_;
  Mat3 current_frame_;
  Mat3 initial_rotation_;
  Mat3 current_rotation_;
  Vec3 pos_;
};

#endif // FRAME_ATTACHMENT_H
