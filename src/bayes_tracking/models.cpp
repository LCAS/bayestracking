/***************************************************************************
 *   Copyright (C) 2006 by Nicola Bellotto                                 *
 *   nbellotto@lincoln.ac.uk                                                    *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "bayes_tracking/models.h"
#include "bayes_tracking/BayesFilter/matSup.hpp"
// #include <angle.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <float.h>


using namespace Bayesian_filter;
using namespace Models;


//*************************************************************************
//                      CV PREDICTION MODEL
//*************************************************************************

CVModel::CVModel(Float wxSD, Float wySD) :
   Linrz_predict_model(x_size, q_size),
   Sampled_predict_model(),
   fx(x_size),
   genn(rnd),
   xp(x_size),
   n(q_size), rootq(q_size),
   m_wxSD(wxSD), m_wySD(wySD)
{
   first_init = true;
   init();
}


void CVModel::init()
{
  Fx.clear();
  // x
  Fx(0,0) = 1.;
  Fx(0,1) = dt;
  Fx(0,2) = 0.;
  Fx(0,3) = 0.;
  // dx
  Fx(1,0) = 0.;
  Fx(1,1) = 1.;
  Fx(1,2) = 0.;
  Fx(1,3) = 0.;
  // y
  Fx(2,0) = 0.;
  Fx(2,1) = 0.;
  Fx(2,2) = 1.;
  Fx(2,3) = dt;
  // dy
  Fx(3,0) = 0.;
  Fx(3,1) = 0.;
  Fx(3,2) = 0.;
  Fx(3,3) = 1.;
  // noise
  q[0] = sqr(m_wxSD); // cov(w_x)
  q[1] = sqr(m_wySD); // cov(w_y)
  G.clear();
  G(0,0) = 0.5*sqr(dt);
  G(1,0) = dt;
  G(2,1) = 0.5*sqr(dt);
  G(3,1) = dt;
}


const FM::Vec& CVModel::f(const FM::Vec& x) const
{  // human model
   fx[0] = x[0] + x[1] * dt;  // x
   fx[1] = x[1];              // dx
   fx[2] = x[2] + x[3] * dt;  // y
   fx[3] = x[3];              // dy
   return fx;
}


void CVModel::update(double dt)
{
  // time interval
  this->dt = dt;

  //     | 0.5dt^2     0    |
  // G = |    dt       0    |
  //     |    0     0.5dt^2 |
  //     |    0        dt   |
  G(0,0) = 0.5*sqr(dt);
  G(1,0) = dt;
  G(2,1) = 0.5*sqr(dt);
  G(3,1) = dt;
}


void CVModel::updateJacobian(const FM::Vec& x) {
   //      | 1  dt 0  0 |
   // Fx = | 0  1  0  0 |
   //      | 0  0  1  dt|
   //      | 0  0  0  1 |
  Fx(0,1) = dt;
  Fx(2,3) = dt;
}


const Vec& CVModel::fw(const FM::Vec& x) const
/*
   * Definition of sampler for additive noise model given state x
   *  Generate Gaussian correlated samples
   * Precond: init_GqG, automatic on first use
   */
{
   if (first_init)
      init_GqG();
                  // Predict state using supplied functional predict model
   xp = f(x);
                  // Additive random noise
   CVModel::genn.normal(n);            // independant zero mean normal
                        // multiply elements by std dev
   for (FM::DenseVec::iterator ni = n.begin(); ni != n.end(); ++ni) {
      *ni *= rootq[ni.index()];
   }
   FM::noalias(xp) += FM::prod(this->G,n);         // add correlated noise
   return xp;
}


void CVModel::init_GqG() const
/* initialise predict given a change to q,G
   *  Implementation: Update rootq
   */
{
   first_init = false;
   for (FM::Vec::const_iterator qi = this->q.begin(); qi != this->q.end(); ++qi) {
      if (*qi < 0)
         error (Numeric_exception("Negative q in init_GqG"));
      rootq[qi.index()] = std::sqrt(*qi);
   }
}



//*************************************************************************
//                   CARTESIAN SUBTRACTION OBSERVATION MODEL
//*************************************************************************

CartesianModel::CartesianModel(Float xSD, Float ySD) :
   Linrz_correlated_observe_model(x_size, z_size),
   Likelihood_observe_model(z_size),
   z_pred(z_size),
   li(z_size)
{
  // Hx = | 1  0  0  0 |
  //      | 0  0  1  0 |
  Hx.clear();
  Hx(0,0) = 1.;
  Hx(1,2) = 1.;
  // noise
  Z.clear();
  Z(0,0) = sqr(xSD);
  Z(1,1) = sqr(ySD);
}

Bayes_base::Float
 CartesianModel::Likelihood_correlated::L(const Correlated_additive_observe_model& model, const FM::Vec& z, const FM::Vec& zp) const
/*
 * Definition of likelihood given an additive Gaussian observation model:
 *  p(z|x) = exp(-0.5*(z-h(x))'*inv(Z)*(z-h(x))) / sqrt(2pi^nz*det(Z));
 *  L(x) the the Likelihood L(x) doesn't depend on / sqrt(2pi^nz) for constant z size
 * Precond: Observation Information: z,Z_inv,detZterm
 */
{
   if (!zset)
      Bayes_base::error (Logic_exception ("BGSubModel used without Lz set"));
               // Normalised innovation
   zInnov = z;
   model.normalise (zInnov, zp);
   FM::noalias(zInnov) -= zp;

   Float logL = scaled_vector_square(zInnov, Z_inv);
   using namespace std;
   return exp(Float(-0.5)*(logL + z.size()*log(2*M_PI) + logdetZ));   // normalized likelihood
}


void CartesianModel::Likelihood_correlated::Lz (const Correlated_additive_observe_model& model)
/* Set the observation zz and Z about which to evaluate the Likelihood function
 * Postcond: Observation Information: z,Z_inv,detZterm
 */
{
   zset = true;
                  // Compute inverse of Z and its reciprocal condition number
   Float detZ;
   Float rcond = FM::UdUinversePD (Z_inv, detZ, model.Z);
   model.rclimit.check_PD(rcond, "Z not PD in observe");
                  // detZ > 0 as Z PD
   using namespace std;
   logdetZ = log(detZ);
}


Bayes_base::Float
 CartesianModel::Likelihood_correlated::scaled_vector_square(const FM::Vec& v, const FM::SymMatrix& V)
/*
 * Compute covariance scaled square inner product of a Vector: v'*V*v
 */
{
   return FM::inner_prod(v, FM::prod(V,v));
}


const FM::Vec& CartesianModel::h(const FM::Vec& x) const
{
  z_pred[0] = x[0];
  z_pred[1] = x[2];
  return z_pred;
};
   

void CartesianModel::updateJacobian(const FM::Vec& x) {
  // nothing to do
}


void CartesianModel::normalise(FM::Vec& z_denorm, const FM::Vec& z_from) const {
}
