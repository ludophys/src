#include <n4-all.hh>
#include <G4ParticleGun.hh>
#include "G4EmStandardPhysics.hh"
#include <G4GenericMessenger.hh>
#define _USE_MATH_DEFINES
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <G4PrimaryParticle.hh>
#include <G4Material.hh>
#include "G4NistManager.hh"
#include <G4String.hh>
#include <G4SystemOfUnits.hh>   // physical units such as `m` for metre
#include <G4Event.hh>           // needed to inject primary particles into an event
#include <G4Box.hh>             // for creating shapes in the geometry
#include <G4Sphere.hh>          // for creating shapes in the geometry
#include <FTFP_BERT.hh>         // our choice of physics list
#include <G4RandomDirection.hh> // for launching particles in random directions

#include <Randomize.hh>

std::ofstream output;
std::ifstream inFile;

std::vector<double> flux;
std::vector<double> energy;

G4String prim_ene[1000];

std::map<std::string, G4double> tempInitialEnergies;
std::map<std::string, G4double> tempInitialTimes;
std::map<std::string, G4String> tempParticleNames;
std::map<std::string, G4String> tempInitialVolumes;

unsigned n_event;
G4int manualTrackID;
G4double MeV_to_amu = 1/931.5;
G4double pi = M_PI;


std::map<G4int, G4double> secondaryInitialEnergies;
std::map<G4int, G4String> secondaryParticleNames;

#include <G4ThreeVector.hh>
#include <cstdlib>

struct my {
  
  G4double       msg_pressure{10 * bar};

 
 
};

auto my_generator(const my& my) {
  return [&](G4Event* event) {
    /*auto neutron = n4::find_particle("neutron");
    auto vertex = new G4PrimaryVertex();
    vertex -> SetPosition(0, 0, -70/2 * cm);
    auto direction = G4ThreeVector{0, 0, 10 * MeV}; // random unit vector which hits scintillator
    auto p = 10.0 * MeV * direction;
    vertex -> SetPrimary(new G4PrimaryParticle(
                           neutron,
                           p.x(), p.y(), p.z()
                         ));*/
    /*auto r = my.particle_dir.mag2() > 0 ? my.particle_dir : G4RandomDirection();
    vertex -> SetPrimary(new G4PrimaryParticle(
                           particle_type,
                           r.x(), r.y(), r.z()
                         ));*/
    //event  -> AddPrimaryVertex(vertex);
    //G4double rand_angle = G4UniformRand();

    //To generate isotropic positions
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> angle_theta(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<> angle_u(-1.0, 1.0);

    G4double theta = angle_theta(gen);   
    G4double u = angle_u(gen);           
    G4double phi = std::acos(u); 
    
    G4double rad = 2 * m;
    //G4double phi = G4UniformRand() * pi;
    //G4double theta = G4UniformRand() * pi;
    //std::cout << "rand is : " << 2 * G4UniformRand() << std::endl;
    G4double x = rad * sin(phi) * cos(theta);
    G4double y = rad * sin(phi) * sin(theta);
    G4double z = rad * cos(phi);  

    //To generate primary neutrons with energies according to the flux
    if (!inFile.is_open()) {
    G4cerr << "Error : impossible to open the flux file" << G4endl;
    return;
    
}
std::string line;
while (std::getline(inFile, line)) {
    std::istringstream iss(line);
    double energy_value;
    double flux_value;
    
    if (!(iss >> flux_value >> energy_value)) {
        G4cerr << "Error at the line : " << line << G4endl;
        continue;
    }
    //G4cout <<"flux value is :" << flux_value <<  "energy value is :" << energy_value << G4endl;
    energy.emplace_back(energy_value);
    flux.emplace_back(flux_value);
}

std::vector<G4double> cdf;
    G4double somme = 0.0;
    for (G4double p : flux) {
        somme += p;
        cdf.push_back(somme);
    }

    
    G4double r = G4UniformRand() * somme;

   
    auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
    size_t index = std::distance(cdf.begin(), it);

    //std::cout <<"n_event is : " << n_event << "energy selected is : " << energy[index] << std::endl;
    
    auto neutron = n4::find_particle("neutron");
    auto particle_gun = std::make_shared<G4ParticleGun>(neutron, 1); // 1 particle
    G4double mass = neutron -> GetPDGMass();
    particle_gun -> SetParticleEnergy(energy[index]*MeV);
    //std::cout << energy[index]*MeV <<" " <<  flux[index] << std::endl;
    particle_gun -> SetParticlePosition({x, y, z});
    particle_gun -> SetParticleMomentumDirection({-x, -y ,-z});
    particle_gun -> GeneratePrimaryVertex(event);


  };
}


n4::sensitive_detector* sensitive_detector(unsigned& n_event) {
  auto process_hits = [&] (G4Step* step) {
    auto pre           = step -> GetPreStepPoint();
    auto post           = step -> GetPostStepPoint();
    auto track = step -> GetTrack();

    auto name = track -> GetDefinition() -> GetParticleName();
    auto trk_id = track -> GetTrackID();
    auto parentID = track -> GetParentID();
    auto step_number = track -> GetCurrentStepNumber();
    auto globalTime = track -> GetGlobalTime();

    auto initial_k_energy = pre  -> GetKineticEnergy();
    
    auto pos = pre -> GetPosition();
    auto pos_initial = track -> GetVertexPosition(); // gets where it was converted, if done so
    auto pos_interaction = post -> GetPosition();
    
    auto cell = pre -> GetTouchable() -> GetVolume() -> GetName();
    auto process = post -> GetProcessDefinedStep() -> GetProcessName();
    
    const G4VProcess* preProcess = pre->GetProcessDefinedStep();
    std::string preProcessName = (preProcess) ? preProcess->GetProcessName() : "LackOfProcess";

    auto k_energy_left = post -> GetKineticEnergy();

    G4double m_i = step->GetTrack()->GetDefinition()->GetPDGMass();

    G4double v_i = pre -> GetVelocity();
    G4double v_f = post -> GetVelocity();

    G4ThreeVector p_i = pre -> GetMomentum();
    G4ThreeVector p_f = post -> GetMomentum();
    G4double E_i = pre->GetTotalEnergy(); // gets kinetic and rest mass energy
    G4double E_f = post->GetTotalEnergy(); // E^2 = p^2 + m^2
    
    // // Compute scattering angle θ
    G4double cosTheta = p_i.unit().dot(p_f.unit());

    G4double M_classical = m_i / (std::pow((v_i / v_f) * cosTheta - 1, 2) - 1); // non-relativistic approach, only really works for protons

    G4ThreeVector p_T = p_i - p_f; // Recoil nucleus momentum , because p_i = p_f + p_T
    G4double pT2 = p_T.mag2(); // Squared magnitude of recoil momentum

    // Compute the estimated target mass using relativistic energy conservation
    // E_i + M = E_f + sqrt(pT^2 + M^2)
    G4double numerator = std::pow((E_i - E_f), 2) - pT2;
    G4double denominator = 2 * (E_i - E_f);

    G4double M = numerator / denominator * MeV_to_amu * -1; // relativistic approach, works for all mass target particles
    int target = std::round(M);

    if (process == "hadElastic" && name == "neutron") {
        output << n_event << " " << trk_id << " " << parentID << " " << name << " " << process << " " << cell << " " << initial_k_energy << " " << k_energy_left << " " << initial_k_energy-k_energy_left << " " << target << " " << globalTime << " " << pos_interaction.x() << " " << pos_interaction.y() << " " << pos_interaction.z() << "\n";
    }

    auto secondaries = step -> GetSecondaryInCurrentStep();
    if (secondaries -> size() > 0) {
        for (int i=0; i<secondaries -> size(); i++) {
            auto secondary = secondaries -> at(i);
            auto particle = secondary -> GetDefinition() -> GetParticleName();
            auto energy_sec = secondary -> GetKineticEnergy();
            auto pos2 = secondary -> GetPosition();
            auto trk_id_sec = secondary -> GetTrackID();
            auto sec_parentID = secondary -> GetParentID();
            
            auto sec_globaltime = secondary -> GetGlobalTime();
            auto vertex = secondary -> GetVertexPosition();
            
            manualTrackID++;
            if (process == "neutronInelastic" || process == "neutronCapture" || name == "gamma") {
                std::string identifier = std::to_string(sec_parentID) + "_" + std::to_string(manualTrackID); // + std::to_string(step_number) + "_" + std::to_string(i);

                tempInitialEnergies[identifier] = energy_sec;
                tempParticleNames[identifier] = particle;
                tempInitialTimes[identifier] = globalTime;
                tempInitialVolumes[identifier] = cell;
            }
        }
    }

    // ignoring neutrons and gammas because each always produces a secondary particle, at least in theory (not in practice)...
    std::string identifier = std::to_string(parentID) + "_" + std::to_string(trk_id);
    if ((tempInitialEnergies.find(identifier) != tempInitialEnergies.end() && track -> GetTrackStatus() == fStopAndKill && name != "neutron" && name != "gamma") || (tempInitialEnergies.find(identifier) != tempInitialEnergies.end() && process == "Transportation" && name != "neutron" && name != "gamma") ) {
        auto sec_energy = tempInitialEnergies[identifier];
        auto sec_particle = tempParticleNames[identifier];
        auto sec_time = tempInitialTimes[identifier];
        auto sec_cell = tempInitialVolumes[identifier];

        auto sec_final_energy = track -> GetKineticEnergy();
        G4String sec_process; // = pre -> GetProcessDefinedStep() -> GetProcessName();

        if (process == "Transportation" || preProcessName == "LackOfProcess" ) {
            sec_process = process;
        } else {
            sec_process = pre -> GetProcessDefinedStep() -> GetProcessName();
        }

        if (sec_energy != sec_final_energy) {
            output << n_event << " " << trk_id << " " << parentID << " " << name << " " << sec_process << " " << sec_cell << " " << sec_energy << " " << sec_final_energy << " " << sec_energy-sec_final_energy << " " << "0" << " " << sec_time << " " << pos_interaction.x() << " " << pos_interaction.y() << " " << pos_interaction.z() << "\n";
        }

        if (process == "Transportation") {
            tempInitialVolumes[identifier] = cell;
            tempInitialEnergies[identifier] = sec_final_energy;
            tempInitialTimes[identifier] = globalTime;
        }
        else {
            tempInitialEnergies.erase(identifier);
            tempParticleNames.erase(identifier);
            tempInitialTimes.erase(identifier);
            tempInitialVolumes.erase(identifier);
        }
    }


    return true; // See https://jacg.github.io/nain4/explanation/process-hits-return-value.html
  };
  return new n4::sensitive_detector{"Detector", process_hits};
}


n4::actions* create_actions(my& my, unsigned& n_event) {
  auto my_event_action = [&] (const G4Event*) {
    
    n_event++;
    
     manualTrackID = 1;
    tempInitialEnergies.clear();
    tempParticleNames.clear();
    tempInitialTimes.clear();
    tempInitialVolumes.clear();
     
     
  };

  return (new n4::        actions{my_generator(my)  })
 -> set( (new n4::   event_action{                  }) -> end(my_event_action) );
}

auto my_geometry(const my& my) {

  //Definition of gas properties
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
  //G4double density, temperature, pressure;
  G4Material* gasMaterial;
  //G4Material* xenon = man->FindOrBuildMaterial("G4_Xe");

  //temperature = 293.15*kelvin; 
  //pressure = 14*100000*pascal;    
  
  auto vessel_length = 1110 * mm;
  auto vessel_diam_ext = 665.2 * mm;
  auto vessel_diam_int = 640 * mm;

  auto closure_vessel_diam_ext = 820 * mm;
  auto closure_vessel_length = 100 * mm;
  auto closure_vessel_diam_int = 560 * mm;
  //auto closure_vessel_z = vessel_length/2 + closure_vessel_length/2;
  auto closure_vessel_z = 475 * mm;
  //std::cout << "closure vessel z is : " << closure_vessel_z << std::endl;
  auto disc_cath_diam = 600.4 * mm;
  auto disc_cath_length = 4 * mm;
  auto disc_cath_z = 0 * mm;

  auto cont_EL_diam = 601.4 * mm;
  auto cont_EL_length = 10.5 * mm;
  auto cont_EL_z = disc_cath_z + 12.5 * mm + cont_EL_length/2;

  auto EL_diam = 525.4 * mm;
  auto EL_length = 0.5 * mm;
  auto EL_z = cont_EL_z - 4.75 * mm;

  auto tube_teflon_diam_int = 570.618 * mm;
  auto tube_teflon_diam_ext = tube_teflon_diam_int + 2 * 3.453 * mm;
  auto tube_teflon_length = 400 * mm;
  auto tube_teflon_z = cont_EL_z + tube_teflon_length/2 + 5 * mm;

  auto gasdrift_diam =  tube_teflon_diam_int;
  auto gasdrift_length = 2 * tube_teflon_length + 4 * 5 * mm;
  auto gasdrift_z = disc_cath_z;

  auto ring_diam_int = 582 * mm;
  auto ring_diam_ext = 588 * mm;
  auto ring_length = 7.01 * mm;
  auto ring_z = disc_cath_z + 418.595 * mm + ring_length/2;

  auto world_rad = gasdrift_diam/2 + 150 * cm;

  auto water  = n4::material("G4_WATER");
  auto air    = n4::material("G4_AIR");
  
  auto vaccum = n4::material("G4_Galactic");

  auto steel  = n4::material("G4_STAINLESS-STEEL");
  auto copper = n4::material("G4_Cu");
  auto teflon = n4::material("G4_TEFLON");
 
  G4double M_Xe = 131.29*g/mole;
  G4Element* Xe = new G4Element("Xenon","Xe", 54., M_Xe);

  //G4double pressure = 5*bar; 
  auto pressure = my.msg_pressure;
  G4double temperature = 293.*kelvin; 
  G4double R =  8.3145*joule/(mole*kelvin);
  G4double density = (pressure * M_Xe) / (R * temperature);
  
  G4Material* xenon = new G4Material("Xenon", density, 1, kStateGas, temperature, pressure); 
  xenon->AddElement(Xe, 1);

  std::cout <<  " density is :" << xenon->GetDensity() / (g/cm3) << " g/cm3" << " Temperature is : " << xenon->GetTemperature()/kelvin << " Kelvin" << " Pressure is : " << (xenon->GetPressure()  / bar)  << " bar" << std::endl;
  
//xenon->GetDensity()
  //temperature = 293.15*kelvin; 
  //pressure = 14*100000*pascal; 
  //gasMaterial = new G4Material("XenonGas", density = xenon->GetDensity(), 1, kStateGas, temperature, pressure);

  auto world  = n4::sphere("world").r(world_rad).volume(vaccum);

  auto vessel = n4::tubs  ("vessel" ).r_inner(vessel_diam_int/2).r_delta(vessel_diam_ext/2 - vessel_diam_int/2).z(vessel_length).place(steel).at(0, 0, disc_cath_z).in(world).now();
  auto closure_vessel1 = n4::tubs  ("closure_vessel1" ).r(closure_vessel_diam_ext/2).z(closure_vessel_length).place(steel).at(0, 0, closure_vessel_z).in(world).now();
  auto closure_vessel2 = n4::tubs  ("closure_vessel2" ).r(closure_vessel_diam_ext/2).z(closure_vessel_length).place(steel).at(0, 0, - closure_vessel_z).in(world).now();


  auto disc_cathode = n4::tubs  ("disc_cathode" ).r(disc_cath_diam/2).z(disc_cath_length).place(steel).at(0, 0, disc_cath_z).in(world).now();
  
  auto EL1 = n4::tubs  ("EL1" ).r(EL_diam/2).z(EL_length).place(steel).at(0, 0, EL_z).in(world).now();
  auto EL2 = n4::tubs  ("EL2" ).r(EL_diam/2).z(EL_length).place(steel).at(0, 0, - EL_z).in(world).now();

  auto cont_EL1 = n4::tubs  ("cont_EL1" ).r_inner(EL_diam/2).r_delta(cont_EL_diam/2 - EL_diam/2).z(cont_EL_length).place(steel).at(0, 0, cont_EL_z).in(world).now();
  auto cont_EL2 = n4::tubs  ("cont_EL2" ).r_inner(EL_diam/2).r_delta(cont_EL_diam/2 - EL_diam/2).z(cont_EL_length).place(steel).at(0, 0, - cont_EL_z).in(world).now();

  auto tube_teflon1 = n4::tubs  ("tube_teflon1" ).r_inner(tube_teflon_diam_int/2).r_delta(tube_teflon_diam_ext/2 - tube_teflon_diam_int/2).z(tube_teflon_length).place(teflon).at(0, 0, tube_teflon_z).in(world).now();
  auto tube_teflon2 = n4::tubs  ("tube_teflon2" ).r_inner(tube_teflon_diam_int/2).r_delta(tube_teflon_diam_ext/2 - tube_teflon_diam_int/2).z(tube_teflon_length).place(teflon).at(0, 0,  - tube_teflon_z).in(world).now();
  
  for (G4int i = 0; i < 27; i++) {
  //for (G4int i = 0; i < 1; i++) {

  G4double ring_delta = 7.99 * mm + ring_length;
  G4double ring_space = ring_delta * i;
  
  auto ring1 = n4::tubs  ("ring1" ).r_inner(ring_diam_int/2).r_delta(ring_diam_ext/2 - ring_diam_int/2).z(ring_length).place(copper).at(0, 0, ring_z - ring_space).in(world).now();
  auto ring2 = n4::tubs  ("ring2" ).r_inner(ring_diam_int/2).r_delta(ring_diam_ext/2 - ring_diam_int/2).z(ring_length).place(copper).at(0, 0, -(ring_z - ring_space)).in(world).now();

}
  auto gasdrift = n4::tubs  ("gasdrift" ).r(gasdrift_diam/2).z(gasdrift_length).sensitive(sensitive_detector(n_event)).place(xenon).in(world).now();

  return n4::place(world).now();
}

int main(int argc, char* argv[]) {
 
  n_event = 0;
  output.open ("output.txt");
  inFile.open("flux/bckg_flux_1500_target_v2_norm.txt");
  my my;

  G4int physics_verbosity = 0;

  // The trailing slash after '/my_geometry' is CRUCIAL: without it, the
  // messenger violates the principle of least surprise.
  auto messenger = new G4GenericMessenger{nullptr, "/my/", "docs: bla bla bla"};
  messenger -> DeclarePropertyWithUnit("msg_pressure"      , "bar"  , my.msg_pressure  );
  /*messenger -> DeclarePropertyWithUnit("bubble_radius"     , "m"  , my.bubble_radius  );
  messenger -> DeclarePropertyWithUnit("socket_rot"        , "deg", my.socket_rot     );
  messenger -> DeclarePropertyWithUnit("particle_energy"   , "keV", my.particle_energy);
  messenger -> DeclareProperty        ("particle"          ,        my.particle_name  );
  messenger -> DeclareProperty        ("particle_direction",        my.particle_dir   );
  messenger -> DeclareProperty        ("physics_verbosity" ,        physics_verbosity );*/

  n4::run_manager::create()
    .ui("neutron_bg", argc, argv)
    .macro_path("macs")
    .apply_early_macro("early-hard-wired.mac")
  
    .physics<FTFP_BERT>(physics_verbosity)
    .geometry([&] { return my_geometry(my); })
    .actions(create_actions(my, n_event))
    
    .run();
    std::cout << "Number of events is : " << n_event << std::endl;


  // Important! physics list has to be set before the generator!

}
