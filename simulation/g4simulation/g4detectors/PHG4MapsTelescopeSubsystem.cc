#include "PHG4MapsTelescopeSubsystem.h"
#include "PHG4MapsTelescopeDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4MapsTelescopeSteppingAction.h"
#include "PHG4Parameters.h"

#include <Geant4/G4GDMLParser.hh>

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4MapsTelescopeSubsystem::PHG4MapsTelescopeSubsystem( const std::string &name, const int lyr, int in_stave_type ):
  PHG4DetectorSubsystem( name, lyr ),
  detector_( 0 ),
  steppingAction_( NULL ),
  eventAction_(NULL),
  layer(lyr),
  stave_type(in_stave_type),
  detector_type(name),
  superdetector(name)
{

  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());

  InitializeParameters();

}

//_______________________________________________________________________
int PHG4MapsTelescopeSubsystem::InitRunSubsystem( PHCompositeNode* topNode )
{
  if (verbosity > 0)
    cout << "PHG4MapsTelescopeSubsystem::Init started" << endl;

  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  // These values are set from the calling macro using the setters defined in the .h file
  if (Verbosity())
  {
    cout << "    create MAPS detector for layer " << layer << endl;
  }
  detector_ = new PHG4MapsTelescopeDetector(topNode, GetParams(), Name());
  detector_->Verbosity(Verbosity());
  detector_->SuperDetector(superdetector);
  detector_->Detector(detector_type);
  detector_->OverlapCheck(CheckOverlap());
  if (Verbosity())
  {
    cout << "    ------ created detector for " << layer << " name " << Name()
         << endl;

    GetParams()->Print();
  }

  if (GetParams()->get_int_param("active"))
  {
    ostringstream nodename;
    if (superdetector != "NONE")
    {
      nodename <<  "G4HIT_" << superdetector;
    }
    else
    {
      nodename <<  "G4HIT_" << detector_type << "_" << layer;
    }
    // create hit list
    PHG4HitContainer* block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
    if ( !block_hits )
    {

      dstNode->addNode( new PHIODataNode<PHObject>( block_hits = new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject" ));

      if (Verbosity())
        cout << PHWHERE << "creating hits node " << nodename.str() << endl;
    }

    PHG4EventActionClearZeroEdep *eventaction = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
    if (GetParams()->get_int_param("absorberactive"))
    {
      nodename.str("");
      if (superdetector != "NONE")
      {
        nodename <<  "G4HIT_ABSORBER_" << superdetector;
      }
      else
      {
        nodename <<  "G4HIT_ABSORBER_" << detector_type << "_" << layer;
      }
      block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
      if ( !block_hits )
      {

        dstNode->addNode( new PHIODataNode<PHObject>( block_hits = new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject" ));

      }
      eventaction->AddNode(nodename.str());
    }
    eventAction_ = dynamic_cast<PHG4EventAction *> (eventaction);
    // create stepping action
    steppingAction_ = new PHG4MapsTelescopeSteppingAction(detector_);
  }
  if (GetParams()->get_int_param("blackhole") && !(GetParams()->get_int_param("active")))
  {
    steppingAction_ = new PHG4MapsTelescopeSteppingAction(detector_);
  }
  return 0;

}

//_______________________________________________________________________
int
PHG4MapsTelescopeSubsystem::process_event( PHCompositeNode * topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    steppingAction_->SetInterfacePointers( topNode );
  }
  return 0;
}


//_______________________________________________________________________
PHG4Detector* PHG4MapsTelescopeSubsystem::GetDetector( void ) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4MapsTelescopeSubsystem::GetSteppingAction( void ) const
{
  return steppingAction_;
}

void
PHG4MapsTelescopeSubsystem::SetDefaultParameters()
{
//  set_default_double_param("place_in_x",0);
//  set_default_double_param("place_in_y",0);
//  set_default_double_param("place_in_z",0);
//  set_default_double_param("rot_in_x",0);
//  set_default_double_param("rot_in_y",0);
//  set_default_double_param("rot_in_z",0);
  set_default_int_param("layer", layer);
  set_default_int_param("stave_type", stave_type);
  set_default_int_param("N_staves", -1);

  set_default_double_param("layer_nominal_radius", NAN);

  set_default_double_param("pixel_x", NAN);
  set_default_double_param("pixel_z", NAN);
  set_default_double_param("pixel_thickness", NAN);
  set_default_double_param("phitilt", NAN);

  set_default_string_param("material", "G4_AIR"); // default - almost nothing
//  set_default_int_param("active", 1);
//  set_default_int_param("absorberactive", 0);
//  set_default_int_param("blackhole", 0);

  set_default_string_param("stave_geometry_file", "ITS.gdml"); // default - almost nothing
}
