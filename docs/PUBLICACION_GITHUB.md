# Publicar Jupiter RT Python en GitHub sin sustituir el README existente

La documentación está pensada para añadirse al repositorio donde **ya existe un README de presentación**. Este paquete no aporta ni sobrescribe `README.md`.

## 1. Contenido recomendado

```text
Jupiter_RT_Python/
├── README.md                  # el que ya mantiene el autor
├── .gitignore
├── pyproject.toml
├── jupiter_rt/
│   ├── ...
│   └── data/temperature_table.npy
└── docs/
    ├── GUIA_TECNICA.md
    ├── FORMATOS.md
    ├── VALIDACION.md
    └── PUBLICACION_GITHUB.md
```

No se incluyen datos de observación, tablas `.opa`, perfiles privados, logs de ejecución, resultados científicos pesados ni entornos Conda. Son archivos diferentes del código fuente que deben gestionarse externamente y cuyo derecho de distribución debe verificarse.

## 2. Añadir los documentos al repositorio existente

Descargar el paquete de **documentación** y descomprimirlo en la raíz del repositorio, o usar el paquete completo de código + docs v0.4.3 si también se desea copiar la versión del código aquí documentada. Si se trabaja sobre un repositorio que ya contiene modificaciones propias, **comparar primero los cambios antes de sobrescribir archivos de código**.

Colocar los cuatro `.md` en `docs/` y `.gitignore` en la raíz.

En el README de presentación existente se puede agregar una línea como:

```markdown
## Documentación técnica

- [Guía técnica y uso](docs/GUIA_TECNICA.md)
- [Archivos de entrada y salida](docs/FORMATOS.md)
- [Validación y estado del proyecto](docs/VALIDACION.md)
```

## 3. Comprobar qué se publicará

En la raíz del repositorio:

```bash
git status
find docs -maxdepth 1 -type f -name '*.md' -print
git check-ignore -v archivo.opa  # comprueba una coincidencia de .gitignore, si procede
```

**No ejecutar `git add .` a ciegas** si no se ha inspeccionado `git status`. Antes de subir, comprobar que no haya rutas personales, credenciales, resultados observacionales no distribuibles ni opacidades enormes.

## 4. Subir los cambios si el repositorio Git ya existe

```bash
cd /ruta/a/Jupiter_RT_Python

git status
git add docs .gitignore pyproject.toml jupiter_rt
git status
git commit -m "Document Jupiter RT v0.4.3 and validation status"
git push
```

El comando anterior **no añade ni modifica el README**, de modo que seguirá estando como estaba. Si también se ha editado el README para incluir enlaces a la guía, incluirlo explícitamente en el `git add`.

Si es un repositorio nuevo todavía sin remoto:

```bash
cd /ruta/a/Jupiter_RT_Python
git init
git branch -M main
git add .gitignore pyproject.toml jupiter_rt docs README.md
git commit -m "Initial public release of Jupiter RT Python v0.4.3"
git remote add origin https://github.com/USUARIO/NOMBRE_REPOSITORIO.git
git push -u origin main
```

Sustituir `USUARIO` y `NOMBRE_REPOSITORIO` por sus valores reales. No ejecutar `git init` ni `git remote add origin` si el repositorio ya está configurado; comprobar antes con `git remote -v`.

## 5. Seguridad y reproducibilidad

- **Opacidades**: se leen desde `--opacity-dir` y no deben subirse al repositorio salvo que exista autorización de distribución.
- **Datos científicos y perfiles**: el `.gitignore` propuesto los excluye por defecto; si hay casos públicos reproducibles, pueden distribuirse de forma deliberada y documentada en otro repositorio, un release o un depósito de datos.
- **Tabla de temperaturas del código**: `jupiter_rt/data/temperature_table.npy` **sí debe versionarse**. El `.gitignore` incluido no excluye los `.npy` del paquete.
- **README**: mantener la descripción existente y enlazar desde allí las nuevas guías.
- **Licencia**: no se ha añadido ninguna automáticamente. Antes de hacer público el código, elegir una licencia acorde con los derechos sobre las traducciones y el material Fortran original; comprobar también posibles condiciones de las tablas de opacidad.
- **Citación**: si se busca uso académico externo, añadir un `CITATION.cff` en el futuro con los autores y referencias correctos. No se han inventado nombres, afiliaciones ni DOIs.

## 6. Publicación sin integración GitHub

Este paquete deja los archivos **preparados para GitHub**, pero no crea un repositorio remoto ni publica nada por sí solo. La transferencia se hace con la interfaz de GitHub o usando los comandos `git` anteriores desde el entorno del autor.
