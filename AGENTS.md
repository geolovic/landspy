# Instrucciones de trabajo para agentes

Estas instrucciones se aplican a todo el repositorio. Las instrucciones
explícitas del usuario tienen prioridad. Comunica los resultados en español
cuando el usuario trabaje en español.

## Estructura y entorno

- El paquete está en `src/landspy`; los tests en `tests`; los scripts de
  medición y sus informes en `benchmarks`.
- Consulta `requirements.txt` y `setup.py` antes de cambiar dependencias.
  GDAL necesita sus bibliotecas nativas y los datos de PROJ. Activa el entorno
  completo antes de ejecutar pruebas, en lugar de usar solo su ejecutable Python.
- Instala el paquete en modo editable con `python -m pip install -e .` cuando
  el entorno todavía no tenga configurada la importación del código local.
- En el entorno cloud usado para este trabajo, si siguen existiendo estas rutas:

  ```bash
  source /workspace/.landspy-conda/etc/profile.d/conda.sh
  conda activate /workspace/.landspy-env
  export MPLBACKEND=Agg
  export MPLCONFIGDIR=/workspace/.landspy-setup/matplotlib
  export OPENBLAS_NUM_THREADS=2
  ```

  En otro entorno, adapta las rutas. No añadas estos directorios locales al
  repositorio ni los conviertas en requisitos de la librería.

## Contratos de DEM y Flow

- `_ix`, `_ixc` y `_zx` son estructuras fundamentales de Flow. Conserva su
  significado, correspondencia entre posiciones y tipos. Comprueba también
  `_nodata_pos` cuando cambies la creación del objeto.
- Los pesos y las distancias de Dijkstra se conservan en float64. Reducirlos
  a float32 puede crear empates y cambiar receptores; no lo hagas como una
  optimización transparente. El experimento documentado cambió tres receptores
  de `tunez` con la configuración por defecto.
- `Flow` usa por defecto `filled=False`, `raw_z=False`, `auxtopo=False`.
  Rellena una sola vez y reutiliza esas elevaciones. Con `filled=True`, no
  debe volver a rellenar. `raw_z=True` conserva las alturas del DEM de entrada
  para `_zx`.
- `DEM.fill()` usa Priority-Flood compilado en `_priority_flood.py`. Por
  defecto devuelve una copia. `inplace=True` sustituye el array del DEM tras
  el cálculo, pero sigue necesitando una copia de trabajo temporal.
- El rellenado actual procesa los valores como alturas, incluido el centinela
  NoData; no restaura una máscara NoData. No cambies esta decisión ni el
  filtrado existente de Flow incidentalmente al optimizar memoria.
- `get_weights()` usa el Dijkstra compilado de `_dijkstra.py`: ocho vecinos,
  longitud 1 o `sqrt(2)`, coste `length * 0.5 * (old_cost + new_cost)` y coste
  inicial cero en las semillas. Mantén el orden de estas operaciones; evita
  `fastmath` si necesitas conservar la equivalencia numérica.
- Fuera de los planos, `np.inf` bloquea la propagación. Los pesos devueltos
  fuera de los planos son `-99999`. Con semillas, los planos inalcanzables
  mantienen distancia infinita. Sin presills se conserva el comportamiento
  de devolución de los costes auxiliares.
- El heap de Dijkstra tiene una entrada por celda en la frontera y permite
  actualizar su prioridad. Sus índices y posiciones usan int32 cuando caben
  e int64 para rásteres mayores. No elimines la alternativa int64.
- MCP ya no se usa en el código de producción de Flow. Los tests todavía usan
  `skimage.graph.MCP_Geometric` como referencia; considera esos tests antes de
  eliminar `scikit-image` de las dependencias.

## Tests y validación

Los tests generan archivos en `tests/data/out` y algunos dependen de los
resultados de los anteriores. Ejecútalos sobre una copia temporal para no
modificar los datos del repositorio. Usa el patrón `test_*.py`: `tests/test.py`
ejecuta una suite y borra resultados al importarse, por lo que no debe entrar
en el descubrimiento automático.

Desde la raíz del repositorio, con el entorno activado:

```bash
python - <<'PY'
from pathlib import Path
import os
import shutil
import tempfile
import unittest

repo = Path.cwd()
with tempfile.TemporaryDirectory(prefix='landspy-tests-') as temporary:
    target = Path(temporary) / 'tests'
    shutil.copytree(repo / 'tests', target)
    (target / 'data' / 'out').mkdir(parents=True, exist_ok=True)
    try:
        os.chdir(target)
        suite = unittest.defaultTestLoader.discover('.', pattern='test_*.py')
        result = unittest.TextTestRunner(verbosity=1).run(suite)
    finally:
        os.chdir(repo)
    raise SystemExit(not result.wasSuccessful())
PY
```

Para cambios en los algoritmos de Flow:

- Ejecuta los tests de la librería. En la implementación documentada pasan
  115 tests; el número puede crecer con nuevos casos.
- Compara la versión anterior con la nueva en `small25`, `tunez` y `jebja30`
  para las ocho combinaciones de `filled`, `raw_z` y `auxtopo` (24 casos).
- Comprueba igualdad de arrays y tipos. Si cambia la ordenación, alinea por
  celda emisora antes de contar cambios de receptores; una diferencia de
  posiciones no equivale necesariamente a una diferencia de drenaje.
- Para el solver, cubre barreras, semillas múltiples y repetidas, zonas
  inalcanzables, rásteres estrechos, crecimiento del heap e índices int64.
- No rebajes tolerancias ni cambies resultados esperados para ocultar una
  regresión. Documenta cualquier diferencia y su efecto en el flujo.

## Medición de memoria y velocidad

Consulta `benchmarks/README.md` para resultados, revisiones y límites de las
comparaciones existentes. No confundas medidas de `fill()` con la creación
completa de Flow ni un ráster aleatorio con un DEM real.

```bash
# Comparación de creación de Flow sobre un DEM sintético pequeño.
python benchmarks/benchmark_flow_memory.py --size 256 --baseline 09e481d

# Perfil del código actual sobre un DEM real.
python benchmarks/profile_flow_memory.py --dem /ruta/DEM.tif --output actual.json

# Referencia MCP sobre exactamente el mismo DEM.
python benchmarks/profile_flow_memory.py --dem /ruta/DEM.tif --revision 09e481d --output mcp.json
```

- El perfil por fases requiere Linux y `/proc`. El muestreo externo cada
  50 ms evita que las extensiones que retienen el GIL bloqueen al muestreador.
- Ejecuta las versiones secuencialmente en procesos separados, con el mismo
  archivo, opciones y entorno. Calienta los kernels antes de medir.
- El tiempo de construcción excluye lectura del DEM, compilación inicial,
  hashing y guardado de arrays. El máximo RSS del sistema corresponde al
  proceso completo; los máximos de cada fase son muestreados y aproximados.
- Los picos de las fases no se suman. Informa el pico global y la fase que
  ahora lo produce; reducir una fase puede desplazar el cuello de botella.
- Usa `--save-arrays /ruta/temporal` si necesitas comparar los arrays completos
  después de medir. Prevé espacio de disco para ambas versiones.
- Conserva los JSON de resultados, dimensiones, dtype, opciones, revisiones
  y comprobaciones. Un único ensayo por versión no establece significación
  estadística de las diferencias de tiempo.
- El tamaño comprimido del TIFF no es el tamaño del DEM en RAM. Usa
  `array.nbytes`; informa también el número de celdas y NoData.
- El sintético histórico tiene 14.480 × 14.480 celdas float32 (799,83 MiB),
  alturas enteras aleatorias de 0 a 1999 y semilla 38. El DEM real suministrado
  tiene 16.627 × 14.448 celdas int16 (458,20 MiB). No los intercambies al
  presentar resultados.
- Mantén los DEM aportados por el usuario, sus vistas y los arrays completos
  fuera del repositorio. Guarda en Git los scripts y resultados agregados,
  como se ha hecho en los benchmarks existentes.

## Cambios y entrega

- Revisa el estado de Git y conserva cambios existentes que no pertenezcan
  a la tarea. No añadas imágenes, datos o archivos temporales por accidente.
- Antes de entregar, ejecuta `git diff --check` y las comprobaciones adecuadas
  al cambio. Describe qué cambió, cómo se verificó y qué falta medir.
- Trabaja y publica en la rama autorizada por el usuario. No mezcles cambios
  en `master` ni crees flujos de aprobación adicionales para acciones ya
  autorizadas. La rama de esta serie de optimizaciones es
  `optimize-priority-flood-fill`; verifica el contexto antes de reutilizarla
  para otra tarea.
