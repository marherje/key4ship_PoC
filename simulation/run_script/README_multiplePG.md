# Multiple Particle Gun simulations

Flujo para generar eventos donde **cada evento contiene varias partículas
primarias simultáneas desde el mismo vértice**, sin ficheros HepMC/HEPEVT de
entrada: las partículas se generan en tiempo de ejecución con el Particle Gun
de DDG4/ddsim.

## Cómo funciona

- `multiplePG_base.py` contiene toda la lógica: registra un
  `Geant4ParticleGun` por partícula a través de
  `SIM.inputConfig.userInputPlugin`. ddsim asigna una interaction *Mask*
  distinta a cada gun y fusiona todas las primarias en **el mismo evento
  Geant4** (`Geant4InteractionMerger` + `Geant4PrimaryHandler`).
- Los steering files de `steer/` son configuración pura (partículas, energías,
  ángulo, vértice, nº de eventos, salida) — no duplican código.
- El ángulo de apertura entre partículas es configurable (`angleDeg`); la
  primera partícula va según +z y las demás se abren ese ángulo (para N>2 se
  reparten en acimut). El vértice es común y configurable.
- Preparado para GENIE u otro generador externo: basta con que
  `configure()` registre un plugin lector (p. ej. `Geant4InputAction` +
  EDM4hep/HepMC) en lugar de los guns; condor y análisis no cambian.
- Por defecto ddsim solo guarda las primarias en `MCParticles`; con
  `keepAllParticles = True` en el dict de configuración se conserva la lista
  completa de secundarias (ficheros mucho más grandes, pero necesario para que
  la comparación primarias vs estado final del análisis sea no trivial).

Casos configurados por defecto (en `launch_multiplePG.sh`):

| Caso | Partículas | Energías |
|---|---|---|
| 1 | μ⁺ (PDG −13) + π⁺ (PDG 211) | 10 + 5 GeV |
| 2 | e⁺ (PDG −11) + π⁺ (PDG 211) | 10 + 5 GeV |

Para añadir un par nuevo: una línea más en el array `pairs` de
`launch_multiplePG.sh` (formato `"particula1:E1:particula2:E2"`).

## Parametrización de la geometría

La geometría se controla desde `simulation/geometry/parameters.yaml` (número
de capas y posición z de cada detector; `auto` = colocación derivada, los
detectores se apilan solos). Para aplicarla:

```bash
cd simulation/geometry
python3 config.py            # parchea SND_compact.xml (backup en .bak)
python3 config.py --dry-run  # ver los cambios sin escribir
python3 config.py --show     # resumen de la geometría actual
```

`config.py` valida antes de escribir (solapes entre detectores, gun aguas
arriba) y `launch_multiplePG.sh` lo invoca automáticamente al inicio, así que
las muestras siempre se producen con la geometría declarada en el YAML. Como
simulación y reconstrucción leen el mismo XML, no hay que tocar ningún
steering ni job de Gaudi. Tras cambiar la geometría hay que re-simular (los
datos existentes dejan de ser consistentes).

## Estructura de salida

```
simulation/run_script/
 ├── steer/   steering files + wrappers de condor generados
 ├── log/     logs separados por configuración (.ddsim.log, outfile_*, errors_*, *.log)
 └── data/    output_<label>.edm4hep.root  (label = <physlist>_SND_<p1>_<E1>GeV_<p2>_<E2>GeV_angle_<a>)
```

## Cómo lanzar

**Una simulación individual desde terminal** (tras generar el steering con
`launch_multiplePG.sh`, o escribiéndolo a mano con la plantilla del docstring
de `multiplePG_base.py`):

```bash
source ../../init_key4hep.sh    # o el entorno de build.sh
ddsim --steeringFile steer/QGSP_BERT_SND_mu+_10GeV_pi+_5GeV_angle_1.0.py
```

**Un conjunto de simulaciones con HTCondor:**

```bash
cd simulation/run_script
./launch_multiplePG.sh                 # genera steering files y hace condor_submit
# un solo job ya generado:
./generic_condor_multiplePG.sh QGSP_BERT_SND_mu+_10GeV_pi+_5GeV_angle_1.0.py \
                               QGSP_BERT_SND_mu+_10GeV_pi+_5GeV_angle_1.0
# sin condor (ejecución local, p. ej. para pruebas):
RUN_LOCAL=1 ./launch_multiplePG.sh
```

**Producción completa con análisis Gaudi:**

```bash
cd simulation/run_script && ./launch_multiplePG.sh     # esperar a los jobs
cd ../../gaudi_jobs/multi_pg_pipeline && ./multi_pg_pipeline.sh
```

## Pipeline de análisis (`gaudi_jobs/multi_pg_pipeline/`)

`job1_mcanalysis.py` corre `MCParticleAnalyzer` (en `gaudi_source/`) sobre un
fichero EDM4hep y escribe `histos_<label>.root` con:

- energía de las primarias (`h_E_primary`), PDG (`h_pdg_primary`)
- distribuciones angulares (`h_theta_primary`, `h_phi_primary`) y ángulo de
  apertura entre las dos primeras primarias (`h_openingAngle`)
- posición del vértice (`h_vtx_x/y/z`)
- nº de MCParticles/primarias/estado-final por evento y comparación
  primarias vs finales (`h_nPrim_vs_nFinal`, `h_Esum_primary`, `h_Esum_final`)

Uso individual:

```bash
INPUT_FILE=../../simulation/run_script/data/output_<label>.edm4hep.root k4run job1_mcanalysis.py
```

Los estudios de calorímetro (energía depositada, hits, shower shape, PID) se
añaden extendiendo `MCParticleAnalyzer.cpp`: las colecciones de SimCalorimeterHit
ya están en los ficheros de entrada.
