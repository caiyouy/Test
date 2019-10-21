// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <vector>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

using namespace Dune;

template<class GridView, class Mapper>
double evolve (const GridView& gridView,
               const Mapper& mapper,
               std :: vector<double>& c,
               const std::function<FieldVector<typename GridView::ctype,GridView::dimensionworld>
                                  (FieldVector<typename GridView::ctype,GridView::dimension>,double)> v,
               const std::function<double
                                   (FieldVector<typename GridView::ctype,GridView::dimension>,double)> b,
               double t)
{
  const int dim = GridView::dimension;
  // allocate a temporary vector for the update
  std::vector<double> update(c.size());
  std::fill(update.begin(), update.end(), 0.0);
  // initialize dt very large
  double dt = std::numeric_limits<double>::max();
  for(const auto& element:elements(gridView))
  {
    auto geometry = element.geometry(); 
    double elementVolume = geometry.volume();
    auto i = mapper.index(element);
    // Will store the total element outflow
    double sumfactor = 0.0;
    for(const auto& intersection:intersections(gridView,element))
    {
      // get geometry type of face
      auto igeo=intersection.geometry();
      auto faceCenter=igeo.center();
      auto velocity=v(faceCenter,t);
      //get the center of the intersecition in local coordinates
      const auto& intersectionReferenceElement=ReferenceElements<double,dim-1>::general
        (intersection.type());
      auto intersectionCenter=intersectionReferenceElement.position(0,0);
      // get normal vector scaled with volume
      auto integrationOuterNormal=intersection.integrationOuterNormal(intersectionCenter);
      // Compute factor occuring in flux fomula
      double factor=velocity*integrationOuterNormal;
      // Upwind Flux
      // Outflow contribution
      update[i]-=c[i]*std::max(0.0,factor)/elementVolume;
      // Inflow contribution
      if(factor<=0) //inflow
      {
        if(intersection.neighbor()) //interior face
        {
          auto j=mapper.index(intersection.outside());
          update[i]-=c[j]*factor/elementVolume;
        }
        if(intersection.boundary())
          update[i]-=b(faceCenter,t)*factor/elementVolume;
      }
      sumfactor+=std::max(0.0,factor);
    }
    dt=std::min(dt,elementVolume/sumfactor);
  }
  // scale dt with safety factor
  dt *=0.99;
  // update the concentratioon vector
  for(int i=0;i<c.size();i++)
    c[i]+=dt*update[i];
  return dt;
}
      
int main(int argc, char** argv)
{
  try{
    // Maybe initialize MPI
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    const int dim=2;
    typedef YaspGrid<dim> GridType;
    GridType grid({1.0,1.0},// upper right corner, the lower left is (0,0)
                  {80,80}); // number of elements per direction
    typedef GridType::LeafGridView GridView;
    GridView gridView=grid.leafGridView();
    MultipleCodimMultipleGeomTypeMapper<GridView,MCMGElementLayout> mapper(gridView);
    std::vector<double> c(mapper.size());
    auto c0=[](auto x){return (x.two_norm()>0.125&& x.two_norm()<0.5)?1.0:0.0;};
    // Initialization
    for (const auto& element: elements(gridView))
    {
      auto geometry=element.geometry();
      auto global=geometry.center();
      c[mapper.index(element)]=c0(global);
    }
    // Write data
    auto vtkWriter=std::make_shared<Dune::VTKWriter<GridView>> (gridView);
    VTKSequenceWriter<GridView> vtkSequenceWriter(vtkWriter,"concentration");
    // Write the initial values
    vtkWriter->addCellData(c,"concentration");
    vtkSequenceWriter.write(0.0);

    double t=0;
    double t_end=0.6;
    int k=0;
    const double saveInterval=0.1;
    double saveStep=0.1;
    auto b=[](const Dune::FieldVector<GridView::ctype,GridView::dimensionworld>& x, double t)
    {
      return 0;
    };
    auto v=[](const Dune::FieldVector<GridView::ctype,GridView::dimensionworld>& x, double t)
    {
      return Dune::FieldVector<double,GridView::dimensionworld> (1.0);
    };
    while (t<t_end)
    {
      double dt=evolve(gridView,mapper,c,v,b,t);
      t+=dt;
      ++k;
      if(t>=saveStep)
      {
        vtkSequenceWriter.write(t);
        saveStep += saveInterval;
      }
      std::cout << "k=" << k<< " t=" << t<< " dt=" << dt <<std::endl;
    }
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e.what() << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
