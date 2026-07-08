"""Shared configuration for multiple-Particle-Gun ddsim steering files.

All the logic for simulating events with SEVERAL primary particles generated
simultaneously from the same vertex lives here. A steering file only has to
provide a small configuration dict:

    import sys
    sys.path.insert(0, "<abs path to simulation/run_script>")
    from DDSim.DD4hepSimulation import DD4hepSimulation
    from multiplePG_base import configure

    SIM = DD4hepSimulation()
    configure(SIM, dict(
        physicsList    = "QGSP_BERT",
        numberOfEvents = 100,
        particles      = [("mu+", 10.0), ("pi+", 5.0)],   # (name, E_total [GeV])
        angleDeg       = 1.0,          # opening angle w.r.t. the first particle
        phiDeg         = 0.0,          # azimuth of the rotation plane (optional)
        position       = (1.0, 1.0, -1000.0),             # common vertex [mm]
        outputFile     = "/abs/path/data/output_<label>.edm4hep.root",
    ))

Adding a new particle pair therefore only requires a new config (one entry in
launch_multiplePG.sh) — no new scripts.

Implementation: one DDG4 ``Geant4ParticleGun`` per particle is registered
through ``SIM.inputConfig.userInputPlugin``. ddsim assigns a distinct
interaction Mask to each plugin and merges all of them into the SAME Geant4
event (Geant4InteractionMerger + Geant4PrimaryHandler), so every event
contains all the primaries simultaneously — this is NOT N independent
single-particle simulations, and no HepMC/HEPEVT intermediate is involved.

GENIE (or any other external generator) can replace the guns later without
touching the launch/condor machinery: make ``configure`` register a reader
plugin (e.g. ``Geant4InputAction`` with an EDM4hep/HepMC file produced by
GENIE) instead of the gun plugins — everything downstream is identical.
"""

import math
import os


def _direction(index, n_particles, angle_deg, phi_deg):
    """Direction of primary `index` (0-based).

    Particle 0 goes along +z (beam). Every other particle is emitted at the
    configured opening angle w.r.t. +z; for more than two particles they are
    spread uniformly in azimuth starting at phiDeg, so N-particle final
    states work out of the box.
    """
    if index == 0:
        return (0.0, 0.0, 1.0)
    alpha = math.radians(angle_deg)
    n_off = max(1, n_particles - 1)
    phi = math.radians(phi_deg) + 2.0 * math.pi * (index - 1) / n_off
    return (math.sin(alpha) * math.cos(phi),
            math.sin(alpha) * math.sin(phi),
            math.cos(alpha))


def _make_gun_plugin(name, particle, energy_gev, direction, position_mm):
    """Return a ddsim userInputPlugin that creates one Geant4ParticleGun.

    Every gun is non-standalone: ddsim assigns it a unique interaction Mask
    and its primaries are merged with the other guns' into one Geant4 event.
    """
    def _plugin(dd4hepSimulation):  # noqa: ARG001 (ddsim passes itself)
        from DDG4 import GeneratorAction, Kernel
        from g4units import GeV
        gun = GeneratorAction(Kernel(), "Geant4ParticleGun/" + name)
        gun.particle     = particle
        gun.Energy       = energy_gev * GeV   # total energy, including mass
        gun.multiplicity = 1
        gun.position     = position_mm        # common primary vertex [mm]
        gun.direction    = direction
        gun.isotrop      = False
        gun.Standalone   = False              # merged by the ddsim input stage
        return gun
    return _plugin


def configure(SIM, cfg):
    """Configure a DD4hepSimulation instance for a multiple-PG run.

    Required cfg keys: physicsList, numberOfEvents, particles, outputFile.
    Optional: angleDeg (default 1.0), phiDeg (0.0), position ((1,1,-1000) mm),
    compactFile (defaults to the repo geometry), randomSeed.
    """
    run_script_dir = os.path.dirname(os.path.abspath(__file__))
    compact = cfg.get(
        "compactFile",
        os.path.join(run_script_dir, "..", "geometry", "SND_compact.xml"))
    compact = os.path.abspath(compact)
    if not os.path.isfile(compact):
        raise RuntimeError("ERROR: geometry file not found: " + compact)

    SIM.runType        = "batch"
    SIM.numberOfEvents = int(cfg["numberOfEvents"])
    SIM.skipNEvents    = 0
    SIM.compactFile    = compact
    SIM._compactFile   = SIM.compactFile
    SIM.outputFile     = os.path.abspath(cfg["outputFile"])
    SIM.physicsList    = cfg["physicsList"]
    if cfg.get("randomSeed") is not None:
        SIM.random.seed = int(cfg["randomSeed"])

    # By default ddsim's Geant4ParticleHandler only stores the primaries (plus
    # decay products) in MCParticles. Enable this to keep the full secondary
    # list (e.g. shower particles) — needed for a meaningful primaries vs
    # final-state comparison, at the cost of much larger files.
    # The default Geant4TCUserParticleHandler must also be disabled: it drops
    # every particle produced outside the "tracker region", which is empty in
    # this geometry (all SND subdetectors are calorimeter-type), so with it
    # active only the primaries would survive.
    if cfg.get("keepAllParticles", False):
        SIM.part.keepAllParticles = True
        SIM.part.userParticleHandler = ""

    # SND-specific sensitive action for the MTC SciFi (same as single-PG sims).
    SIM.action.mapActions["MTC"] = "SND_SciFiAction"

    # One gun per primary, all from the same vertex, merged into one event.
    SIM.enableGun = False
    particles = list(cfg["particles"])
    if not particles:
        raise RuntimeError("cfg['particles'] must contain at least one particle")
    angle_deg = float(cfg.get("angleDeg", 1.0))
    phi_deg   = float(cfg.get("phiDeg", 0.0))
    position  = tuple(cfg.get("position", (1.0, 1.0, -1000.0)))

    plugins = []
    for i, (pname, energy_gev) in enumerate(particles):
        direction = _direction(i, len(particles), angle_deg, phi_deg)
        plugins.append(_make_gun_plugin(
            "Gun%d_%s" % (i, pname.replace("+", "plus").replace("-", "minus")),
            pname, float(energy_gev), direction, position))
    SIM.inputConfig.userInputPlugin = plugins

    print("[multiplePG] %d primaries/event from vertex %s, opening angle %.3f deg:"
          % (len(particles), position, angle_deg))
    for i, (pname, energy_gev) in enumerate(particles):
        print("[multiplePG]   %d: %-6s E=%.2f GeV dir=%s"
              % (i, pname, energy_gev,
                 tuple(round(c, 6) for c in _direction(i, len(particles),
                                                       angle_deg, phi_deg))))
    return SIM
