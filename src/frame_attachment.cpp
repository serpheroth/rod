#include "frame_attachment.h"
#include "vector_lib.h"


FrameAttachment::FrameAttachment(PBDElasticRod* rod,
                                 Mat3& initial_frame, Mat3& attached_frame, Vec3& pos)
  : rod_(rod)
  , initial_frame_(initial_frame)
  , rod_frame_(attached_frame)
  , pos_(pos)
{
  dj::MulMatrixRightTransposed3x3<real>((real (*)[3]) &rod_frame_[0][0],
                                        (real (*)[3]) &initial_frame_[0][0],
                                        (real (*)[3]) &initial_rotation_[0][0]);
  current_frame_ = initial_frame_;
  dj::MulMatrix3x3<real>((real (*)[3]) &initial_rotation_[0][0],
                        (real (*)[3]) &current_frame_[0][0],
                        (real (*)[3]) &current_rotation_[0][0]);
}

int FrameAttachment::UpdateFrame(FrameAttachment::FrameUpdator updator)
{
  updator(*this);
  dj::MulMatrix3x3<real>((real (*)[3]) &initial_rotation_[0][0],
                        (real (*)[3]) &current_frame_[0][0],
                        (real (*)[3]) &current_rotation_[0][0]);
}

void FrameAttachment::EnforceFrameAttachment(Vec3& p0, Vec3& p1, Vec3& pe, real length)
{
  length *= real(0.5);
  Mat3 d1;
  PBDElasticRod::ComputeOneMaterialFrame(p0, p1, pe, d1);
  // Darboux vector
  Mat3& d0 = current_rotation_;
  real x = 1 + d0[0] * d1[0] + d0[1] * d1[1] + d0[2] * d1[2];
  Vec3 darboux_vector;
  x = real(2.0) / (length * x);
  {
    int i = 0, j = 2, k = 1;
    darboux_vector[i] = d0[j] * d1[k] - d0[k] * d1[j];
  }
  {
    int i = 1, j = 0, k = 2;
    darboux_vector[i] = d0[j] * d1[k] - d0[k] * d1[j];
  }
  {
    int i = 2, j = 1, k = 0;
    darboux_vector[i] = d0[j] * d1[k] - d0[k] * d1[j];
  }
  darboux_vector *= x;
  Mat3 dbjpi[3][3];
  // Frame derivative
  ComputeFrameDerivative(p0, p1, pe, d1,
                         dbjpi[0][0], dbjpi[0][1], dbjpi[0][2],
                         dbjpi[1][0], dbjpi[1][1], dbjpi[1][2],
                         dbjpi[2][0], dbjpi[2][1], dbjpi[2][2]);

  Mat3 constraint_jacobian[3];
  Mat3& omega_pb = constraint_jacobian[0];
  Mat3& omega_pc = constraint_jacobian[1];
  Mat3& omega_pe = constraint_jacobian[2];
  const int permutation[3][3] = {
    0, 2, 1,
    1, 0, 2,
    2, 1, 0
  };
  for (int c = 0; c < 3; ++c) {
    const int i = permutation[c][0];
    const int j = permutation[c][1];
    const int k = permutation[c][2];
    // pb
    {
      Vec3 term1(real(0), real(0), real(0));
      Vec3 term2(real(0), real(0), real(0));
      Vec3 tmp(real(0), real(0), real(0));
      // first term
      dj::MulVecMatrix3x3<real>(d0[k](), (real (*)[3]) &dbjpi[j][0], term1());
      dj::MulVecMatrix3x3<real>(d0[j](), (real (*)[3]) &dbjpi[k][0], tmp());
      term1 = term1 - tmp;
      // second term
      for (int n = 0; n < 3; ++n) {
        dj::MulVecMatrix3x3<real>(d0[n](), (real (*)[3]) &dbjpi[n][0], tmp());
        term2 = term2 + tmp;
      }
      omega_pb[i] = (term1) + (real(0.5) * darboux_vector[i] * length) * (term2);
      omega_pb[i] *= (-x * rod::bending_twising_stiffness[i]);
    }
    // pc
    {
      Vec3 tmp(real(0), real(0), real(0));
      // first term
      Vec3 term1(real(0), real(0), real(0));
      dj::MulVecMatrix3x3<real>(d0[k](), (real (*)[3]) &dbjpi[j][1], term1());
      dj::MulVecMatrix3x3<real>(d0[j](), (real (*)[3]) &dbjpi[k][1], tmp());
      term1 = term1 - tmp;
      // second term
      Vec3 term2(real(0), real(0), real(0));
      for (int n = 0; n < 3; ++n) {
        dj::MulVecMatrix3x3<real>(d0[n](), (real (*)[3]) &dbjpi[n][0], tmp());
        term2 = term2 + tmp;
      }
      omega_pc[i] = (term1) + (real(0.5) * darboux_vector[i] * length) * (term2);
      omega_pc[i] *= (-x * rod::bending_twising_stiffness[i]);
    }
    // pe
    {
      Vec3 tmp(real(0), real(0), real(0));
      // first term
      Vec3 term1(real(0), real(0), real(0));
      dj::MulVecMatrix3x3<real>(d0[k](), (real (*)[3]) &dbjpi[j][2], term1());
      dj::MulVecMatrix3x3<real>(d0[j](), (real (*)[3]) &dbjpi[k][2], tmp());
      term1 -= tmp;
      // second term
      Vec3 term2(real(0), real(0), real(0));
      for (int n = 0; n < 3; ++n) {
        dj::MulVecMatrix3x3<real>(d0[n](), (real (*)[3]) &dbjpi[n][2][0][0], tmp());
        term2 += tmp;
      }
      omega_pe[i] = (term1) + (real(0.5) * darboux_vector[i] * length) * (term2);
      omega_pe[i] *= (-x * rod::bending_twising_stiffness[i]);
    }
  }


  Vec3 constraint_value(rod::bending_twising_stiffness[0] * darboux_vector[0],
                        rod::bending_twising_stiffness[1] * darboux_vector[1],
                        rod::bending_twising_stiffness[2] * darboux_vector[2]);
  Vec3* verts[] = {&p0, &p1, &pe};
  for (int c = 0; c < 3; ++c) {
    real factor = 0;
    for (int i = 0; i < 3; ++i) {
      factor += constraint_jacobian[i][c] * constraint_jacobian[i][c];
    }
    factor = -constraint_value[c] / factor;
//        factor *= 1e-1f;
    for (int i = 0; i < 3; ++i) {
      *(verts[i]) += (factor * constraint_jacobian[i][c]);
    }
  }
}

void FrameAttachment::EnforcePosAttachment(Vec3& p)
{
  p = pos_;
}

inline void FrameAttachment::ComputeFrameDerivative(Vec3 &p0, Vec3 &p1, Vec3 &p2, Mat3& d,
                                                    Mat3& d1p0, Mat3& d1p1, Mat3& d1p2,
                                                    Mat3& d2p0, Mat3& d2p1, Mat3& d2p2,
                                                    Mat3& d3p0, Mat3& d3p1, Mat3& d3p2)
{

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // d3pi
  Vec3 p01 = p1 - p0;
  real length_p01 = p01.Magnitude();
  d3p0[0] = d[2][0] * d[2];
  d3p0[1] = d[2][1] * d[2];
  d3p0[2] = d[2][2] * d[2];
  d3p0[0][0] -= 1;
  d3p0[1][1] -= 1;
  d3p0[2][2] -= 1;
  d3p0[0] *= (1 / length_p01);
  d3p0[1] *= (1 / length_p01);
  d3p0[2] *= (1 / length_p01);

  d3p1[0] = real(-1.0) * d3p0[0];
  d3p1[1] = real(-1.0) * d3p0[1];
  d3p1[2] = real(-1.0) * d3p0[2];

  d3p2[0].Fill(real(0.0));
  d3p2[0].Fill(real(0.0));
  d3p2[0].Fill(real(0.0));
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // d2pi
  Vec3 p02 = p2 - p0;
  Vec3 p01_cross_p02;
  dj::Cross3(p01(), p02(), p01_cross_p02());
  real length_cross = p01_cross_p02.Magnitude();
  Mat3 mat(d[1][0] * d[1],
           d[1][1] * d[1],
           d[1][2] * d[1]);
  mat[0][0] -= 1;
  mat[1][1] -= 1;
  mat[2][2] -= 1;
  mat[0] *= (real(-1.0) / length_cross);
  mat[1] *= (real(-1.0) / length_cross);
  mat[2] *= (real(-1.0) / length_cross);
  //  Vec3 p2p1 = p2 - p1;
  Mat3 product_matrix;
  GetCrossProductMatrix(p2 - p1, product_matrix);
  dj::MulMatrix3x3<real>(&mat[0][0], &product_matrix[0][0], &d2p0[0][0]);
  GetCrossProductMatrix(p0 - p2, product_matrix);
  dj::MulMatrix3x3<real>(&mat[0][0], &product_matrix[0][0], &d2p1[0][0]);
  GetCrossProductMatrix(p1 - p0, product_matrix);
  dj::MulMatrix3x3<real>(&mat[0][0], &product_matrix[0][0], &d2p2[0][0]);
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // d1pi
  Mat3 product_mat_d3 = GetCrossProductMatrix(d[2]);
  Mat3 product_mat_d2 = GetCrossProductMatrix(d[1]);
  Mat3 m1, m2;

  dj::MulMatrix3x3<real>(&product_mat_d3[0][0], &d2p0[0][0], &m1[0][0]);
  dj::MulMatrix3x3<real>(&product_mat_d2[0][0], &d3p0[0][0], &m2[0][0]);
  d1p0 = m2 - m1;

  dj::MulMatrix3x3<real>(&product_mat_d3[0][0], &d2p1[0][0], &m1[0][0]);
  dj::MulMatrix3x3<real>(&product_mat_d2[0][0], &d3p1[0][0], &m2[0][0]);
  d1p1 = m2 - m1;

  dj::MulMatrix3x3<real>(&product_mat_d3[0][0], &d2p2[0][0], &d1p2[0][0]);
  d1p2[0] *= real(-1);
  d1p2[1] *= real(-1);
  d1p2[2] *= real(-1);
}

