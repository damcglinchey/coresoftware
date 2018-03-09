#include "MvtxStandaloneTracking.h"

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrDefUtil.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TDecompSVD.h>
#include <TMath.h>

#include <iostream>
#include <utility>
#include <stdio.h>

MvtxStandaloneTracking::MvtxStandaloneTracking()
  : window_x_(10)
  , window_z_(10)
  , verbosity_(0)
{

}

MvtxStandaloneTracking::~MvtxStandaloneTracking()
{
}

void
MvtxStandaloneTracking::RunTracking(PHCompositeNode* topNode, MvtxTrackList &trklst)
{

  trklst.clear();

  //------
  //--- get cluster container
  //------
  clusters_ = findNode::getClass<TrkrClusterContainer>(topNode, "TrkrClusterContainer");
  if (!clusters_)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TrkrClusterContainer" << std::endl;
    return;
  }


  //------
  //--- associate clusters
  //------
  AssociateClusters(trklst);

  if ( verbosity_ > 0 )
    std::cout << PHWHERE << " Finished associating clusters. N candidates:" << trklst.size() << std::endl;
  //------
  //--- fit tracks
  //------
  for ( unsigned int itrk = 0; itrk < trklst.size(); itrk++)
  {
    TrackFitXY(trklst.at(itrk));
    TrackFitZY(trklst.at(itrk));
  } // itrk

  if ( verbosity_ > 1 )
    PrintTrackCandidates(trklst);


  //------
  // --- choose best tracks
  //------

  // --- done
  return;
}

void
MvtxStandaloneTracking::AssociateClusters(MvtxTrackList &trklst)
{

  // --- get clusters
  TrkrDefUtil util;

  std::multimap<int, TrkrCluster*> lyrclustermap;
  TrkrClusterContainer::ConstRange clusrange = clusters_->GetClusters();
  for ( TrkrClusterContainer::ConstIterator iter = clusrange.first;
        iter != clusrange.second;
        ++iter)
  {
    TrkrDefs::cluskey ckey = (iter->second)->GetClusKey();
    int lyr = util.GetLayer(ckey);

    lyrclustermap.insert(std::make_pair(lyr, iter->second));
  }

  // --- make track candidates

  // layer 0 clusters
  for ( auto iter0 = lyrclustermap.lower_bound(0);
        iter0 != lyrclustermap.upper_bound(0);
        ++iter0)
  {
    double mx = (iter0->second)->GetX();
    double mz = (iter0->second)->GetZ();
    // layer 1 clusters
    for ( auto iter1 = lyrclustermap.lower_bound(1);
          iter1 != lyrclustermap.upper_bound(1);
          ++iter1)
    {
      mx += (iter1->second)->GetX();
      mz += (iter1->second)->GetZ();

      // layer 2 clusters
      for ( auto iter2 = lyrclustermap.lower_bound(2);
            iter2 != lyrclustermap.upper_bound(2);
            ++iter2)
      {
        mx += (iter2->second)->GetX();
        mz += (iter2->second)->GetZ();

        for ( auto iter3 = lyrclustermap.lower_bound(3);
              iter3 != lyrclustermap.upper_bound(3);
              ++iter3)
        {
          mx += (iter3->second)->GetX();
          mz += (iter3->second)->GetZ();

          // average
          mx /= 4.;
          mz /= 4.;

          // check windows
          if ( fabs((iter0->second)->GetX() - mx) > window_x_ ||
               fabs((iter0->second)->GetZ() - mz) > window_z_ )
            continue;

          if ( fabs((iter1->second)->GetX() - mx) > window_x_ ||
               fabs((iter1->second)->GetZ() - mz) > window_z_ )
            continue;

          if ( fabs((iter2->second)->GetX() - mx) > window_x_ ||
               fabs((iter2->second)->GetZ() - mz) > window_z_ )
            continue;

          if ( fabs((iter3->second)->GetX() - mx) > window_x_ ||
               fabs((iter3->second)->GetZ() - mz) > window_z_ )
            continue;

          // make candidate
          MvtxTrack trk;
          (trk.ClusterList).push_back(iter0->second);
          (trk.ClusterList).push_back(iter1->second);
          (trk.ClusterList).push_back(iter2->second);
          (trk.ClusterList).push_back(iter3->second);

          trklst.push_back(trk);
        } // layer 3
      } // layer 2
    } // layer 1
  } // layer 0

}


TVectorD
MvtxStandaloneTracking::SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L)
{
  // Simple generalized linear least-squares fitter.
  // Solve y(X) = beta'*X + error(L) for beta (vector of regression coefs.)
  // L is the inverse covariance matrix for y.
  // Least-squares solution involves solving X' L X beta = X' L y

  TMatrixD XT(X); XT.T();
  TMatrixD A = XT * L * X;
  TVectorD b = XT * L * y;

  // Now solve A*beta = b using SVD. Decompose A = U S V'.
  TDecompSVD svd(A);
  TMatrixD UT = svd.GetU(); UT.T();

  // Construct Moore-Penrose pseudoinverse of S
  TVectorD s = svd.GetSig();
  TMatrixD Sd(s.GetNrows(), s.GetNrows());
  for (int i = 0; i < s.GetNrows(); i++)
    Sd(i, i) = s(i) > 0 ? 1. / s(i) : 0.;

  TVectorD beta = svd.GetV() * Sd * UT * b;

  return beta;
}

void
MvtxStandaloneTracking::TrackFitXY(MvtxTrack &trk)
{
  // Longitudinal/polar component
  // Fit track using a straight line x(y') = x0 + c*y'.
  // Assigns residuals, z-intercept, and polar angle.

  // m = # measurements; n = # parameters.
  int m = (trk.ClusterList).size(), n = 2;

  TMatrixD X(m, n);
  TMatrixD Cinv(m, m);
  TVectorD y(m);
  TVectorD x(m);

  for (int iclus = 0; iclus < m; iclus++)
  {
    TrkrCluster *clus = (trk.ClusterList).at(iclus);

    x(iclus) = clus->GetX();
    y(iclus) = clus->GetY();

    X(iclus, 0) = 1;
    X(iclus, 1) = y(iclus);

    Cinv(iclus, iclus) = clus->GetSize(0, 0);
  }

  TVectorD beta = SolveGLS(X, x, Cinv);
  double x0 = beta(0), c = beta(1);

  trk.m_xy = c;
  trk.b_xy = x0;

  double chi2  = 0;
  for (int iclus = 0; iclus < m; iclus++)
  {
    TrkrCluster *clus = (trk.ClusterList).at(iclus);

    double cx = clus->GetX();
    double cy = clus->GetY();

    double px = cy * c + x0;

    double dx = cx - px;
    trk.dx.push_back(dx);

    chi2 += pow(dx, 2) / clus->GetError(0, 0);
  }
  chi2 /= double(m);
  trk.chi2_xy = chi2;

  return;
}

void
MvtxStandaloneTracking::TrackFitZY(MvtxTrack &trk)
{
  // Longitudinal/polar component
  // Fit track using a straight line z(y') = z0 + c*y'.
  // Assigns residuals, z-intercept, and polar angle.

  // m = # measurements; n = # parameters.
  int m = (trk.ClusterList).size(), n = 2;

  TMatrixD X(m, n);
  TMatrixD Cinv(m, m);
  TVectorD y(m);
  TVectorD z(m);

  for (int iclus = 0; iclus < m; iclus++)
  {
    TrkrCluster *clus = (trk.ClusterList).at(iclus);

    z(iclus) = clus->GetZ();
    y(iclus) = clus->GetY();

    X(iclus, 0) = 1;
    X(iclus, 1) = y(iclus);

    Cinv(iclus, iclus) = clus->GetSize(2, 2);
  }

  TVectorD beta = SolveGLS(X, z, Cinv);
  double z0 = beta(0), c = beta(1);

  trk.m_zy = c;
  trk.b_zy = z0;

  double chi2 = 0;
  for (int iclus = 0; iclus < m; iclus++)
  {
    TrkrCluster *clus = (trk.ClusterList).at(iclus);

    double cz = clus->GetZ();
    double cy = clus->GetY();

    double pz = cy * c + z0;

    double dz = cz - pz;
    trk.dz.push_back(dz);

    chi2 += pow(dz, 2) / clus->GetError(2, 2);
  }
  chi2 /= double(m);
  trk.chi2_zy = chi2;

  return;
}

void
MvtxStandaloneTracking::PrintTrackCandidates(MvtxTrackList &trklst)
{
  std::cout << "===================================================" << std::endl;
  std::cout << "== " << PHWHERE << " Found " << trklst.size() << " Track Candidates" << std::endl;
  std::cout << "===================================================" << std::endl;
  for ( unsigned int i = 0; i < trklst.size(); i++)
  {
    std::cout << "== " << i << std::endl;
    for ( unsigned int j = 0; j < trklst.at(i).ClusterList.size(); j++)
    {
      std::cout << "    clus " << j
                << " key:0x" << std::hex << trklst.at(i).ClusterList.at(j)->GetClusKey() << std::dec
                << " (" << trklst.at(i).ClusterList.at(j)->GetX()
                << ", " << trklst.at(i).ClusterList.at(j)->GetY()
                << ", " << trklst.at(i).ClusterList.at(j)->GetZ()
                << ")"
                << " dx:" << trklst.at(i).dx.at(j)
                << " dz:" << trklst.at(i).dz.at(j)
                << std::endl;
    }
    std::cout << "    xy"
              << " m:" << trklst.at(i).m_xy
              << " b:" << trklst.at(i).b_xy
              << " chi2:" << trklst.at(i).chi2_xy
              << std::endl;
    std::cout << "    zy"
              << " m:" << trklst.at(i).m_zy
              << " b:" << trklst.at(i).b_zy
              << " chi2:" << trklst.at(i).chi2_zy
              << std::endl;
  } // i
  std::cout << "===================================================" << std::endl;
}

