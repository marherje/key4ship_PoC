# Optimization — barrido de geometrías

Genera variantes de la geometría del SND, simula cada una y la pasa por la
cadena de reconstrucción, y recoge una figura de mérito por variante. Local y
secuencial: una variante detrás de otra, y el fallo de una no aborta las demás.

> **¿Barrido manual u optimización?** `controller.py` es la herramienta para
> ejecutar a mano un conjunto de variantes que tú escribes en un YAML. El bucle
> de **optimización bayesiana multiobjetivo** — que propone las geometrías él
> solo, las lanza en paralelo a HTCondor y aprende de los resultados — vive en
> [`mobo/`](mobo/README.md) y reutiliza exactamente esta misma cadena
> (`config.py` → `run_sim.py` → jobs de `Analysis/`).

La baseline commiteada, `simulation/geometry/SND_compact.xml`, **nunca se
toca**: cada variante es un compact nuevo en
`Simulation/geometry/variants/varNNNNN.xml`.

## Requisito de entorno

```bash
source init_key4ship.sh     # desde la raíz del repo; hace falta ddsim y k4run
```

Sin él, el controller aborta con un mensaje claro. `--dry-run` es la excepción:
sólo renderiza y valida geometría, y con python3 + PyYAML del sistema basta.

## Layout

```
Optimization/
├── controller.py                     # orquestador
├── results.csv                       # una fila por variante ejecutada
├── Simulation/
│   ├── geometry/                     # maquinaria de variantes
│   │   ├── SND_compact_template.xml  # el compact con placeholders XnameX
│   │   ├── parameters_template.yaml  # valores base de esos placeholders
│   │   ├── variants.yaml             # qué variantes generar
│   │   ├── config.py                 # renderiza + valida una variante
│   │   ├── make_variants.py          # genera un set entero sin simular
│   │   └── variants/                 # compacts generados
│   ├── run_scripts/run_sim.py        # ddsim de una variante
│   └── output/variants/varNNNNN/     # steering.py, ddsim.log, edm4hep
└── Analysis/
    ├── job1_overlay.py               # copias parametrizadas de mu_pi_pipeline
    ├── job4_tracking.py
    ├── job5_rntuple.py
    ├── compute_fom.py                # nhits del SiPad -> fom.txt
    ├── run_analysis.sh               # re-run manual de la cadena
    └── variants/varNNNNN/            # events/tracks/ShipHits.root, fom.txt, logs
```

Las definiciones comunes (elements, materials, los XML de detector,
`parse_geometry.py` y la baseline) se quedan en `simulation/geometry` y se
comparten; las variantes las referencian con rutas relativas.

## Uso

```bash
python3 controller.py                      # todas las variantes del spec
python3 controller.py --dry-run            # sólo geometría, no escribe nada
python3 controller.py --only thin_W        # un subconjunto (repetible)
python3 controller.py --events 5           # simulaciones más cortas
python3 controller.py --only thin_W --skip-sim   # reaprovecha la simulación
python3 controller.py --spec mi_barrido.yaml
```

`--skip-sim` no asigna ids nuevos: busca el `varNNNNN` de una ejecución previa
casando el `name` contra el `params.yaml` de cada directorio de salida, así que
la variante tiene que haberse corrido antes.

Para sólo generar compacts, sin simular:

```bash
cd Simulation/geometry && python3 make_variants.py variants.yaml
```

## Formato del spec

```yaml
base: parameters_template.yaml    # opcional, es el default
variants:
  thin_W:      { SiPad_WThickness: 5, SiPad_NLayers: 32 }
  short_sipad: { SiPad_dim_z: 200, SiPad_NLayers: 12, SiTarget_NLayers: 136 }
```

Cada entrada hereda todos los parámetros de la base y sobrescribe sólo lo que
lista. Sobrescribir una clave que no existe en la base es un error (guarda
contra erratas). Los parámetros disponibles y sus límites están documentados en
`Simulation/geometry/parameters_template.yaml`.

Dos avisos que ahorran depuración:

- `SiTarget_dim_z = 1700 mm − SiPad_dim_z`. Al cambiar `SiPad_dim_z` casi
  siempre hay que recontar `SiTarget_NLayers`, o la variante falla por
  `SiTarget_NLayers_max` aunque lo que tocabas fuese el SiPad.
- `SiPad_layer_gap: auto` dimensiona el hueco de aire de cada capa para que las
  `SiPad_NLayers` capas llenen `SiPad_dim_z` exacto y los planos sensibles
  queden equidistantes. Es el default; ponle un número para fijarlo a mano.

## Resultados

`results.csv` lleva una fila por variante: `var_id`, `name`, los parámetros de
entrada, `nhits_SiPad` y `status`. El `status` dice en qué paso falló
(`geometry_error`, `sim_failed`, `job1_failed`, …), y el log correspondiente
está en el directorio de la variante.
