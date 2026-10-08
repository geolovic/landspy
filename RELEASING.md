# Publicación de landspy

## PyPI: configuración inicial del propietario

En https://pypi.org/manage/project/landspy/settings/publishing/, añade un
Trusted Publisher de GitHub con estos campos exactos:

- Owner: `geolovic`
- Repository: `landspy`
- Workflow filename: `publish.yml`
- Environment: `pypi`

En GitHub, crea el entorno `pypi` en Settings → Environments. La publicación
no necesita guardar un token PyPI en GitHub ni compartir credenciales.

Antes de publicar, ejecuta los tests descritos en AGENTS.md y:

```bash
python -m pip install build twine
python -m build
python -m twine check --strict dist/*
```

Después de configurar el Trusted Publisher, publica la etiqueta:

```bash
git tag -a v1.4.1 -m "landspy 1.4.1"
git push origin v1.4.1
```

El workflow construye wheel y sdist, verifica la versión y la exclusión de
benchmarks, y publica ambos archivos en PyPI. No reutilices una versión ya
publicada. `workflow_dispatch` sobre master solo construye; para reintentar
una publicación selecciona la etiqueta `v1.4.1` o reejecuta su workflow.

## conda-forge

La receta v1 está en `conda-recipe/recipe.yaml`. El SHA256 local se valida
durante la preparación, pero antes de enviar la receta debe sustituirse por
el del sdist realmente publicado: GitHub Actions reconstruye el archivo y
puede producir un hash diferente. Después de publicar en PyPI:

```bash
python conda-recipe/update_source.py
```

El script obtiene de PyPI la URL y SHA256 del sdist de la versión declarada
en la receta. No necesita autenticación.

Prueba con `rattler-build build --recipe conda-recipe/recipe.yaml -c conda-forge`.
Para probar antes de publicar, utiliza una copia temporal de la receta y
cambia `source.url` por `source.path` apuntando al sdist local; elimina
`source.sha256` de esa copia.

Para solicitar la incorporación:

1. Crea un fork de https://github.com/conda-forge/staged-recipes.
2. En una rama nueva, copia `recipe.yaml` y `run_test.py` a `recipes/landspy/`.
3. Ejecuta `conda-smithy recipe-lint --conda-forge recipes/landspy`.
4. Abre un PR con esos dos archivos. `geolovic` figura como mantenedor.
5. Espera la revisión y las compilaciones de conda-forge. Su aprobación
   crea el feedstock y publica el paquete; la receta en este repositorio
   por sí sola no lo publica en el canal.

Para futuras versiones, cambia `setup.py`, el contexto `version` de la receta
y CHANGELOG.md. Publica primero en PyPI, actualiza el SHA256, y actualiza el
feedstock existente cuando conda-forge ya lo haya creado.
