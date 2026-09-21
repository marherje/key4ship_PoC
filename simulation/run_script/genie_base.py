"""Shared configuration for GENIE-driven ddsim steering files.

Counterpart of ``multiplePG_base.py`` for external-generator input instead
of particle guns, exactly as anticipated in that module's docstring: a
HepMC file produced from GENIE output replaces the gun plugins, and nothing
downstream (digitization, tracking, RNTuple pipelines) has to change.

The HepMC file itself is produced by ``genie_gst_to_hepmc.py`` from a GENIE
"gst" summary tree (see that script's docstring for the event-building
logic and reference-frame notes). A steering file only has to provide a
small configuration dict:

    import sys
    sys.path.insert(0, "<abs path to simulation/run_script>")
    from DDSim.DD4hepSimulation import DD4hepSimulation
    from genie_base import configure

    SIM = DD4hepSimulation()
    configure(SIM, dict(
        physicsList = "QGSP_BERT",
        hepmcFile   = "/abs/path/data/genie_allflavours_0_500.hepmc",
        outputFile  = "/abs/path/data/output_genie_allflavours_0_500.edm4hep.root",
    ))
"""

import os


def configure(SIM, cfg):
    """Configure a DD4hepSimulation instance to read primaries from a HepMC file.

    Required cfg keys: hepmcFile, outputFile.
    Optional: physicsList (default QGSP_BERT), numberOfEvents (default -1,
    i.e. consume the whole HepMC file -- see DD4hepSimulation.validate(),
    a negative count is only valid together with inputFiles), compactFile
    (defaults to the repo geometry), randomSeed, keepAllParticles.
    """
    run_script_dir = os.path.dirname(os.path.abspath(__file__))
    compact = cfg.get(
        "compactFile",
        os.path.join(run_script_dir, "..", "geometry", "SND_compact.xml"))
    compact = os.path.abspath(compact)
    if not os.path.isfile(compact):
        raise RuntimeError("ERROR: geometry file not found: " + compact)

    hepmc_file = os.path.abspath(cfg["hepmcFile"])
    if not os.path.isfile(hepmc_file):
        raise RuntimeError(
            "ERROR: HepMC file not found: %s "
            "(generate it with genie_gst_to_hepmc.py)" % hepmc_file)

    SIM.runType        = "batch"
    SIM.numberOfEvents = int(cfg.get("numberOfEvents", -1))
    SIM.skipNEvents    = int(cfg.get("skipNEvents", 0))
    SIM.compactFile    = compact
    SIM._compactFile   = SIM.compactFile
    SIM.outputFile     = os.path.abspath(cfg["outputFile"])
    SIM.physicsList    = cfg.get("physicsList", "QGSP_BERT")
    if cfg.get("randomSeed") is not None:
        SIM.random.seed = int(cfg["randomSeed"])

    # See multiplePG_base.configure() for why both flags are needed together.
    if cfg.get("keepAllParticles", False):
        SIM.part.keepAllParticles = True
        SIM.part.userParticleHandler = ""

    # SND-specific sensitive action for the MTC SciFi (same as PG sims).
    SIM.action.mapActions["MTC"] = "SND_SciFiAction"

    # GENIE final-state particles replace the particle gun entirely: ddsim
    # auto-selects Geant4InputAction/HepMC3FileReader from the .hepmc
    # extension (see DD4hepSimulation.py, POSSIBLEINPUTFILES).
    SIM.enableGun = False
    SIM.inputConfig.userInputPlugin = []
    SIM.inputFiles = [hepmc_file]

    # Per-event chatter off unless asked for. At INFO ddsim writes ~1 KB per
    # event (GenerationInit, the HepMC reader, hepmc4, PrimaryHandler,
    # ParticleHandler, Geant4Output2EDM4hep -- seven lines each), which is
    # 10 MB per 10k-event trial and, over a 200-trial campaign, more disk than
    # the results. WARNING (4) keeps the geometry summary and every problem.
    # `verbose: true` in the steering dict restores the default.
    if not cfg.get("verbose", False):
        SIM.printLevel = 5
        SIM.output.inputStage = 4
        SIM.output.kernel = 4
        SIM.output.part = 4
        SIM.output.random = 4
        SIM.output.geometry = 4

    print("[genie_base] reading primaries from %s (numberOfEvents=%d, skip=%d)"
          % (hepmc_file, SIM.numberOfEvents, SIM.skipNEvents))
    return SIM
