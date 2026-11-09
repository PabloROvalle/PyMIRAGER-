# Guía técnica de Jupiter RT Python

**Versión documentada:** 0.4.3  
**Estado del proyecto:** versión modular operativa; validación progresiva frente a implementaciones Fortran originales.  
**Fecha de corte de esta guía:** 9 de noviembre de 2026.

> Esta guía complementa el `README.md` de presentación del repositorio; no lo sustituye. Explica el funcionamiento real de la versión 0.4.3, no funcionalidades planificadas como si ya estuvieran implementadas.

## Índice

1. [Objetivo y alcance](#1-objetivo-y-alcance)
2. [Estructura del código](#2-estructura-del-código)
3. [Requisitos e instalación](#3-requisitos-e-instalación)
4. [Preparación de un caso](#4-preparación-de-un-caso)
5. [Ejecución](#5-ejecución)
6. [Modelo físico y numérico](#6-modelo-físico-y-numérico)
7. [Inversión y regularización](#7-inversión-y-regularización)
8. [Presets disponibles](#8-presets-disponibles)
9. [Jacobianos](#9-jacobianos)
10. [Lectura e interpretación de resultados](#10-lectura-e-interpretación-de-resultados)
11. [Diagnóstico y solución de problemas](#11-diagnóstico-y-solución-de-problemas)
12. [Qué está validado y qué falta](#12-qué-está-validado-y-qué-falta)

Consulta además [formatos de entrada y salida](FORMATOS.md), [validación científica](VALIDACION.md) y [publicación en GitHub](PUBLICACION_GITHUB.md).

## 1. Objetivo y alcance

`jupiter_rt` es una implementación **Python** de un modelo de emisión térmica infrarroja y de inversión de perfiles atmosféricos de Júpiter, derivada de varias versiones Fortran. Permite recuperar temperatura, abundancias moleculares, aerosoles y algunas moléculas especiales, manteniendo configuraciones heredadas y comparando alternativas de jacobianos moleculares.

El intervalo general de trabajo es **500–2000 cm⁻¹** (aproximadamente 5–20 µm). El modelo está construido para **361 niveles de presión**; todavía no se ha generalizado a otras dimensiones verticales. La rejilla interna es de **0,001 cm⁻¹**, en bloques de 10 cm⁻¹. La respuesta instrumental se aplica mediante convolución gaussiana según `resolution.txt`.

La publicación del código **no incluye las tablas binarias de opacidad `.opa`**: se suministran por separado en el entorno de ejecución. Tampoco exige compilación de rutinas Fortran en tiempo de ejecución: los cálculos de absorción inducida por colisiones (CIA) están escritos en Python.

## 2. Estructura del código

```text
Jupiter_RT_Python/
├── jupiter_rt/
│   ├── __init__.py
│   ├── __main__.py      # python -m jupiter_rt
│   ├── cli.py           # argumentos y presets
│   ├── constants.py     # constantes físicas y planetarias
│   ├── geometry.py      # geometría, alturas y gravedad
│   ├── io.py            # lectura de casos y escritura de resultados
│   ├── opacity.py       # opacidades binarias y rejillas espectrales
│   ├── cia_python.py    # componentes CIA H2-H2, H2-He y H2-CH4
│   ├── continuum.py     # profundidad óptica del continuo
│   ├── forward.py       # emisión térmica, convolución y Jacobianos
│   ├── jacobians.py     # alternativas moleculares
│   ├── inverse.py       # iteración y regularización
│   └── data/
│       └── temperature_table.npy  # tabla necesaria, incluida
├── docs/
│   ├── GUIA_TECNICA.md
│   ├── FORMATOS.md
│   ├── VALIDACION.md
│   └── PUBLICACION_GITHUB.md
└── pyproject.toml
```

**No es necesario editar código para cambiar de molécula** si el `.fuel` está correctamente configurado y existe su tabla `.opa` en la ubicación esperada. El `preset` establece valores de regularización, filtro y convolución; **no cambia las especies seleccionadas en el `.fuel`**.

## 3. Requisitos e instalación

- Linux (entorno de validación) con Python **3.10 o superior**.
- Dependencias declaradas: `numpy>=1.24`, `scipy>=1.10`.
- Acceso a los archivos de cada caso y a las tablas `.opa` correspondientes.
- Memoria suficiente para la rejilla espectral y los jacobianos.

Instalación editable en el entorno Conda que se utilizará para ejecutar:

```bash
conda activate jwst_neptune
cd ~/MIRAGER/Jupiter_RT_Python
python3 -m pip install -e .
python3 -c "import jupiter_rt; print(jupiter_rt.__file__)"
python3 -m jupiter_rt --help
```

La ruta `~/MIRAGER/Jupiter_RT_Python` es la del servidor usado durante la validación; **no es un requisito del programa**.

Tras `pip install -e .`, el módulo puede ejecutarse desde la carpeta del caso. El paquete también declara el comando de consola `jupiter-rt`; la documentación utiliza `python3 -m jupiter_rt` para hacer explícito qué intérprete se emplea.

## 4. Preparación de un caso

Un caso contiene los siguientes archivos:

```text
mi_caso/
├── jupiter_ch4.fuel
├── jupiter_ch4.pta
├── jupiter_ch4.spe
├── resolution.txt
└── corrch4.txt
```

Para aerosoles (`cloud 1` o `cloud 2`) debe existir también, respectivamente, `nh3_clouds.txt` o `ch4_hazes.txt`. Las tablas de absorción no tienen que estar en la misma carpeta; se indican con `--opacity-dir`.

Los nombres exactos están definidos por el lector y deben respetarse. Los detalles de columnas, flags y escalados se explican en [FORMATOS.md](FORMATOS.md).

**Convención fundamental del `.pta`:** primera columna presión (bar), segunda temperatura (K), y **una columna de VMR por gas, en el mismo orden que `.fuel`**. No existe columna de aerosoles: `taucloud` comienza como un vector de ceros. Para las moléculas especiales HCN, CO₂, H₂O y C₆H₆, el gas especial va **el último** tanto en `.fuel` como en `.pta`.

La versión 0.4.3 acepta asimismo, solo para casos sin molécula especial, una columna final adicional si **todos** sus valores son exactamente cero; la ignora por compatibilidad con ciertos archivos antiguos.

### Opacidades

El directorio indicado en `--opacity-dir` debe contener un archivo `jupiter_ch4_<gas>.opa` por cada gas listado en `.fuel` (incluidos los fijos). La correspondencia es literal con los nombres de `.fuel`: **`c2h3` se mantiene como identificador de archivo aunque en la interpretación científica del caso se hable de C₂H₂**. No se renombra automáticamente.

Ejemplo de ubicación en el servidor de validación:

```text
/obs/tfouchet/data/opacite/jupiter/jwst
```

No se debe usar como ubicación operativa la ruta heredada `/poubelle/...` que pueda aparecer en logs antiguos de Fortran.

## 5. Ejecución

### Verificación de los archivos (sin opacidades)

```bash
python3 -m jupiter_rt --case . --preset co2 --check-inputs
```

Verifica lectura de `.fuel`, `.pta`, `.spe` y resolución instrumental. **No** lee las tablas `.opa` ni realiza un cálculo radiativo; tampoco certifica por sí solo el contenido de `corrch4.txt`.

### Caso con recuperación de abundancia normal

```bash
python3 -m jupiter_rt \
  --case . \
  --opacity-dir /obs/tfouchet/data/opacite/jupiter/jwst \
  --preset abundance \
  --jacobian original \
  --output resultados_abundance
```

### Caso con molécula especial (CO₂)

```bash
python3 -m jupiter_rt \
  --case . \
  --opacity-dir /obs/tfouchet/data/opacite/jupiter/jwst \
  --preset co2 \
  --jacobian original \
  --output resultados_co2
```

Para C₆H₆ se dispone de `c6h6` y `c6h6_fortran`, con distintas regularizaciones. Para cambiar de especie hay que reemplazar **los datos del caso**, no solamente el argumento `--preset`.

### Solo modelo directo

```bash
python3 -m jupiter_rt \
  --case . \
  --opacity-dir /obs/tfouchet/data/opacite/jupiter/jwst \
  --preset temperature \
  --forward-only \
  --output resultados_forward
```

### Procesos de larga duración

```bash
nohup python3 -u -m jupiter_rt \
  --case . \
  --opacity-dir /obs/tfouchet/data/opacite/jupiter/jwst \
  --preset co2 \
  --output resultados_co2 \
  > log_co2.txt 2>&1 &
```

Consultar `tail -f log_co2.txt`; comprobar el proceso con `ps -u "$USER" -o pid,etime,cmd | grep '[j]upiter_rt'`.

### Parámetros de línea de comandos

| Argumento | Efecto |
|---|---|
| `--case PATH` | Directorio de archivos de entrada; **obligatorio**. |
| `--opacity-dir PATH` | Directorio de tablas `.opa`; si se omite, intenta usar `--case`. |
| `--preset NAME` | Configuración predeterminada; por omisión `auto`. |
| `--jacobian METHOD` | `original` (por defecto), `analytical` o `finite_difference`. |
| `--output PATH` | Directorio de resultados, por defecto `resultados_python`. Se interpreta respecto al directorio actual si es relativo. |
| `--check-inputs` | Revisión de entradas sin calcular opacidades. |
| `--forward-only` | Solo síntesis espectral. |
| `--options FILE.json` | Modifica opciones de inversión. |
| `--opacity-metadata FILE.json` | Define rejillas alternativas de opacidades. |
| `--endian little|big` | Endianness de `.opa`; por defecto `little`. |
| `--max-samples N` | Límite de muestras espectrales; por defecto 160 000. |
| `--without-ch4-legacy-scale` | Desactiva el factor histórico `7/6` aplicado al primer gas. |

### Opciones de inversión personalizadas

`--options` recibe un JSON. Por ejemplo:

```json
{
  "max_iterations": 30,
  "rtol_chi2": 0.001,
  "c_temperature": 0.75,
  "c_special": 100.0,
  "factor": 3.0,
  "stop_mode": "combined"
}
```

El archivo sobrescribe solo las claves indicadas; el resto procede del preset. `c_*` controla la correlación vertical, y `factor` interviene en el escalado de la matriz de covarianza. **No conviene comparar retrievals con regularizaciones diferentes atribuyendo todas las diferencias al modelo directo o al jacobiano.**

## 6. Modelo físico y numérico

Secuencia general implementada:

1. `io.py`: carga perfiles de presión, temperatura, gases y observaciones.
2. `geometry.py`: calcula altura, gravedad dependiente de latitud y geometría de salida.
3. `opacity.py`: lee coeficientes de opacidad binarios, interpolados en temperatura, para cada gas.
4. `cia_python.py` y `continuum.py`: añaden CIA **H₂–H₂, H₂–He y H₂–CH₄**. No se utiliza una biblioteca Fortran externa para la operación estándar.
5. `forward.py`: acumula profundidades ópticas, transmisiones y emisión térmica por capas; convoluciona con la respuesta instrumental y calcula jacobianos cuando procede.
6. `inverse.py`: regulariza, calcula incrementos y actualiza iterativamente temperatura, VMR y/o aerosoles.

La radiancia monocromática térmica se construye usando la función de Planck y diferencias de transmisión entre capas. La temperatura de cada capa es el promedio de las temperaturas de dos niveles vecinos. No se ha reformulado la física del continuo para ganar velocidad si eso alteraba los resultados heredados.

La convolución instrumental es gaussiana y usa una anchura efectiva configurable en cada preset, junto con la resolución de `resolution.txt`.

### Rejillas de absorción

Los valores por defecto de `opacity.py` son:

| Tipo de tabla | Rango de bloques binarios |
|---|---|
| Gases normales | 500–2000 cm⁻¹ |
| HCN | 680–750 cm⁻¹ |
| CO₂ | 640–710 cm⁻¹ |
| H₂O | 500–2000 cm⁻¹ |
| C₆H₆ | 640–710 cm⁻¹ |

Cada bloque tiene 10 cm⁻¹, 10 000 muestras y paso 0,001 cm⁻¹. Los rangos anteriores describen la **cobertura de las tablas binarias**, no necesariamente el intervalo observado de un caso. Si las tablas reales tienen otra estructura, se puede usar `--opacity-metadata` con la especificación apropiada; no hay detección automática de su formato físico.

## 7. Inversión y regularización

El vector de parámetros está construido según los flags de `.fuel`: temperatura cuando existe un gas con flag 1 o 3; VMR de gases con flag 2 o 3; aerosoles cuando `cloud` es distinto de cero.

Se minimiza un desajuste espectral regularizado, trabajando con una correlación gaussiana en logaritmo de presión:

\[
S_{ij}(c)=\exp\!\left[-\frac12\left(\frac{\ln p_i-\ln p_j}{c}\right)^2\right].
\]

El código asigna la longitud de correlación `c_temperature`, `c_molecule`, `c_special` o `c_aerosol` a cada bloque. Calcula un factor adaptativo `alpha` que incorpora la opción `factor`, las incertidumbres observacionales y la norma de sensibilidad del jacobiano. Después resuelve un sistema espectral de covarianza, obtiene las actualizaciones `delta`, las incertidumbres y los averaging kernels.

La **actualización molecular es logarítmica**: `VMR_new = VMR_old * exp(delta)`. Por eso un gas recuperado debe tener VMR estrictamente positiva antes de invertir. La temperatura se actualiza de manera aditiva y la opacidad de aerosoles de manera aditiva con límite inferior cero.

Existen máscaras verticales heredadas basadas en tangentes hiperbólicas, centradas alrededor del índice 90. Su comportamiento exacto depende de `filter_mode`:

- `legacy_abundance` y `legacy_co2`: filtro troposférico para el primer bloque de parámetros y estratosférico para el segundo.
- `legacy_aerosol` y `legacy_c6h6`: filtro troposférico para el primer bloque.
- `legacy_hcn`: multiplica el incremento del primer bloque por `scale_hcn_delta` (0,1 por defecto).
- `legacy_h2o` y `none`: no introducen una máscara adicional en la función `apply_filters` de esta versión.

La anchura de la transición implementada actualmente es de **15 niveles**. Esta decisión se conserva en la versión 0.4.3 por continuidad con las pruebas científicas; no se debe afirmar que todos los filtros reproducen literalmente cada versión Fortran.

### Criterios de parada

- `fortran`: comprueba la variación relativa absoluta del χ² respecto al valor **actual**.
- `combined`: además exige que `max|delta|` disminuya y sea inferior a `tol_step` (0,3 por defecto).
- El máximo de iteraciones depende del preset.

Por pantalla se imprime el **cambio relativo firmado**

\[
(\chi_i^2-\chi_{i-1}^2)/\chi_{i-1}^2,
\]

con `N/A` para la primera iteración. Este número de diagnóstico tiene un denominador distinto del utilizado internamente en el test `rtol_chi2` y no representa necesariamente la convergencia por sí solo.

## 8. Presets disponibles

Valores de `cli.py` y `RetrievalOptions` en v0.4.3:

| Preset | Convolución `gaussian_width` | Longitud de correlación relevante | Filtro | Parada |
|---|---:|---|---|---|
| `temperature` | 1,000 | `c_T=0,75` | `none` | `fortran` |
| `abundance` | 0,925 | `c_T=3`, `c_mol=3` | `legacy_abundance` | `fortran`, máx. 20 |
| `abundance_split` | 0,925 | `c_T=0,75`, `c_mol=3` | `legacy_abundance` | `fortran`, máx. 20 |
| `aerosol` | 1,000 | `c_aerosol=0,75` | `legacy_aerosol` | `combined` |
| `hcn` | 1,000 | `c_special=100` | `legacy_hcn` | `combined` |
| `co2` | 1,200 | `c_special=100` | `legacy_co2` | `combined` |
| `h2o` | 1,000 | `c_special=100` | `legacy_h2o` | `combined` |
| `c6h6` | 0,925 | `c_special=100` | `legacy_c6h6` | `combined` |
| `c6h6_fortran` | 0,925 | `c_T=2`, `c_special=2` | `legacy_c6h6` | `fortran`, máx. 20 |
| `auto` | 0,925 | según tipo de parámetro | auto | `combined` |

Salvo los valores explícitos de la tabla, se heredan los defaults: `c_T=0,75`, `c_mol=3`, `c_special=100`, `c_aerosol=0,75`, `factor=3`, `max_iterations=15`, `rtol_chi2=0,01`. En todos los presets listados el `factor` general es 3.

**Importante:** los nombres `fortran` y `legacy` identifican intenciones de compatibilidad, pero no garantizan bit a bit el mismo resultado que todos los programas Fortran originales. Por ejemplo, el jacobiano térmico y determinados filtros difieren de algunas versiones históricas.

## 9. Jacobianos

La opción `--jacobian` afecta al **jacobiano molecular**, no al térmico ni al de aerosoles.

- `original` (predeterminado): reproduce la fórmula heredada sobre una rejilla espectral submuestreada (`nu[::5]`). Es el valor de referencia para mantener continuidad de resultados.
- `analytical`: evalúa una derivada molecular alternativa a partir de la transmisión, utilizando la rejilla completa.
- `finite_difference`: evalúa diferencias centradas de las **expresiones algebraicas de transmisión**, no dos ejecuciones completas e independientes del modelo directo con atmósfera perturbada. Se utiliza como comprobación interna de la implementación analítica.

El jacobiano térmico distribuye entre niveles adyacentes la contribución de la derivada de Planck por capas. El jacobiano de aerosoles mantiene la formulación original del proyecto, pendiente de una validación independiente completa.

Se ha comprobado que los métodos molecular `analytical` y `finite_difference` prácticamente coinciden entre sí para el caso molecular estudiado. Esto **no** sustituye una comprobación de diferencias finitas del modelo directo completo. Para detalles y cifras, consultar [VALIDACION.md](VALIDACION.md).

## 10. Lectura e interpretación de resultados

Una inversión genera los archivos siguientes en la carpeta de `--output`:

| Archivo | Contenido |
|---|---|
| `jupiter_ch4.res` | Número de onda, espectro observado, sintético y error (`4` columnas). |
| `jupiter_ch4.inv` | Presión, temperatura recuperada, perfiles de gases seleccionados, aerosoles si activos, incertidumbres; consultar cabecera. |
| `jupiter_ch4.sum` | Iteración final, χ² medio y DOFS por parámetro. |
| `jupiter_ch4.avg` | Para cada nivel/columna de la matriz de averaging kernels, **índice 1-based del máximo**; no guarda la matriz 361×361 completa. |
| `diagnostics.json` | χ² final, iteraciones, parámetros, longitudes de correlación, método de jacobiano e historial iterativo. |

`--forward-only` escribe el espectro y perfiles sin los diagnósticos completos de la inversión.

El estadístico mostrado es el **χ² medio**, definido como el promedio de los cuadrados de los residuos normalizados por la incertidumbre; no es automáticamente el «χ² reducido» con un número explícito de grados de libertad residuales.

Los DOFS son trazas de los averaging kernels. Los máximos de `jupiter_ch4.avg` permiten estudiar dónde alcanza el máximo cada columna, pero **no** su anchura, su estructura completa ni el número de niveles verticales independientes por sí solos.

**Atención al escalado:** `io.py` multiplica por los factores declarados en `.fuel`; de manera adicional aplica `7/6` al **primer gas** por compatibilidad histórica, salvo `--without-ch4-legacy-scale`, y **×50 exclusivamente a C₆H₆**. En `.inv`, para los gases recuperados, se divide por el factor `scale` de `.fuel`, **no** por los factores históricos `7/6` y `50`. Por tanto, comprueba la convención antes de comparar una columna `.inv` directamente con los números crudos del `.pta`.

## 11. Diagnóstico y solución de problemas

| Problema | Causa habitual | Comprobación |
|---|---|---|
| `No module named jupiter_rt` | Paquete no instalado en el intérprete activo | `python3 -m pip install -e .` dentro del repositorio y verificar `python3 -c 'import jupiter_rt'`. |
| `FileNotFoundError: jupiter_ch4.fuel` | `--case` apunta a carpeta errónea | `ls /ruta/caso/jupiter_ch4.fuel`; `--case .` significa directorio actual. |
| Se invierte `c2h3` al elegir preset `c6h6` | Archivos de entrada incorrectos | Comprobar `jupiter_ch4.fuel`, los puntos de `.spe` y la línea `Retrieved parameters` de `--check-inputs`. |
| `Missing opacity file` | `.opa` ausente o nombre distinto | Comprobar `<opacity-dir>/jupiter_ch4_<nombre>.opa`. |
| Tamaño de opacidad incorrecto | Tabla con estructura diferente | Verificar formato real y metadatos; no alterar a ciegas los bloques. |
| Columna extra en `.pta` | Confusión entre VMR y aerosol | La última columna de un caso especial es su VMR; no existe columna de aerosoles. |
| Perfil molecular no positivo | Recuperación logarítmica requiere VMR positiva | Revisar los valores iniciales de moléculas con flags 2/3. |
| Cambios grandes de perfiles con χ² similar | Correlaciones/regul. diferentes, baja información vertical | Comparar presets, DOFS, incertidumbres y kernels, no solo χ². |
| Falta `corrch4.txt` | CIA no puede calcularse | Colocar el archivo en `--case` y ejecutar un caso real. |

## 12. Qué está validado y qué falta

- **Continuo CIA:** comparación numérica directa de componentes con el Fortran original.
- **Temperatura:** muy buena concordancia de recuperación Python–Fortran en el caso de referencia estudiado.
- **Abundancia normal C₂H₂ (`c2h3`) con ajuste térmico:** muy buena concordancia con la referencia.
- **C₆H₆:** el caso correcto se ejecuta y el espectro calculado coincide muy bien con Fortran; las discrepancias verticales de abundancia se atribuyen, según la configuración usada, a regularizaciones distintas. Se considera operativamente aceptado para continuar, **sin reclamar equivalencia exacta de perfiles**.
- **CO₂:** siguiente caso planificado para validar la rama compartida de moléculas especiales. **Todavía no validado** al cierre de esta guía.
- **HCN y H₂O:** no comparados individualmente; que CO₂ funcione no verifica automáticamente sus tablas de opacidad, sus filtros o sus espectros particulares.
- **Aerosoles:** se ha conservado el tratamiento heredado; falta una comparación completa de extremo a extremo documentada en esta versión.
- **Jacobianos:** `analytical` y `finite_difference` son internamente consistentes; falta la perturbación de todo el modelo directo y una comparación de matrices completas de averaging kernels.

Los números, condiciones y advertencias de interpretación están en [VALIDACION.md](VALIDACION.md). Las futuras validaciones deberían registrar versión, preset, `.fuel`, rejilla espectral, regularización, jacobiano, número de iteraciones y referencia Fortran, para evitar comparar casos distintos.
