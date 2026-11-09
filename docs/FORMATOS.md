# Formatos de datos, factores de escala y archivos de salida

**Versión:** Jupiter RT Python 0.4.3. Esta es la especificación **implementada por el lector actual**, no una descripción genérica de todos los archivos Fortran históricos.

## 1. Estructura de un caso

```text
case_directory/
  jupiter_ch4.fuel
  jupiter_ch4.pta
  jupiter_ch4.spe
  resolution.txt
  corrch4.txt
  nh3_clouds.txt     # solo si cloud 1
  ch4_hazes.txt      # solo si cloud 2
```

Los archivos `.opa` se almacenan habitualmente en una carpeta independiente y se pasan con `--opacity-dir`.

## 2. `jupiter_ch4.fuel`

El lector elimina comentarios posteriores a `!` o `#` y líneas vacías. Las líneas útiles son, en este orden:

1. Latitud, en grados.
2. Coseno del ángulo de incidencia.
3. Coseno del ángulo de emergencia.
4. `cloud N` con `N` igual a 0, 1 o 2.
5. Una línea por gas con tres campos: **nombre, factor de escala, flag**.

Esquema ilustrativo para un caso especial:

```text
1.003336
0.994964
0.984629
cloud 0
ch5 1.0 1
ch3e 1624.0 0
nh3 1.0 0
ph3 1.0 0
c2h3 1.0 0
c2h6 1.0 0
c6h6 1.0 2
```

Los números del ejemplo únicamente ilustran el orden; utilizar los valores científicos de cada observación. Las etiquetas son nombres de ficheros de opacidad, no nomenclatura química normalizada.

### Flags

| Flag | Tratamiento |
|---:|---|
| 0 | Gas fijo; se usa su perfil y opacidad, no se invierte. |
| 1 | Activa recuperación de **temperatura** (por la presencia de algún gas con flag 1); el gas no se invierte por VMR. |
| 2 | Recupera la **VMR del gas**. |
| 3 | Activa temperatura y recuperación de **VMR del gas**. |

La temperatura aparece como un único bloque de parámetros aunque haya varios gases con flag 1 o 3. No basta con seleccionar un preset para activar una molécula: el flag de `.fuel` debe estar configurado correctamente.

La molécula especial, si existe, **debe ser la última** de la lista de gases.

### Nubes

- `cloud 0`: sin inversión de aerosoles.
- `cloud 1`: activa aerosoles y lee `nh3_clouds.txt`.
- `cloud 2`: activa aerosoles y lee `ch4_hazes.txt`.

En todos los casos `taucloud` empieza en cero. Los archivos de nubes proporcionan la **máscara vertical** (361 filas, presión más máscara), no una columna de opacidad inicial en `.pta`.

## 3. `jupiter_ch4.pta`

Matriz numérica ASCII con **361 filas**:

```text
pressure_bar  temperature_K  VMR_gas_1  VMR_gas_2  ...  VMR_gas_N
```

La columna 1 contiene presión decreciente (bar). La columna 2, temperatura (K). Las columnas posteriores corresponden exactamente a la lista de gases de `.fuel`.

Ejemplo de **orden de columnas**, no de valores:

```text
P | T | ch5 | ch3e | nh3 | ph3 | c2h3 | c2h6 | c6h6
```

**Regla permanente:** NO hay columna de aerosoles. Para HCN, CO₂, H₂O y C₆H₆, **la última columna es su VMR**.

La versión 0.4.3 admite una columna adicional final de ceros **solo en casos sin molécula especial** para leer archivos heredados; nunca se interpreta como opacidad de aerosoles.

### Aplicación de escalas al leer

1. Se cargan las VMR de todas las columnas del `.pta`.
2. Por compatibilidad con el código histórico se multiplica el **primer perfil molecular** por `7/6` (se puede desactivar con `--without-ch4-legacy-scale`).
3. Se multiplica cada perfil por su factor numérico de `.fuel`.
4. Si el gas se llama exactamente `c6h6`, se multiplica **además por 50**, una sola vez.

**HCN, CO₂ y H₂O no reciben el factor ×50.** Tampoco reciben ×5 por ser moléculas especiales.

La regla `7/6` se aplica al primer gas listado, incluso si su etiqueta no es literalmente `ch4`; por tanto se debe mantener el orden de gases esperado por el modelo de referencia.

### Interpretación de las VMR en `.inv`

El escritor divide los gases recuperados por el factor declarado en `.fuel`, pero **no deshace** los factores adicionales `7/6` o `50`. No comparar sin revisar el escalado la última columna del `.pta` con la de `.inv`.

## 4. `jupiter_ch4.spe`

Archivo ASCII cuyo primer registro contiene el número `N` de puntos observados, seguido de **N filas** con:

```text
wavenumber_cm-1   radiance_observed   sigma_radiance
```

`N` debe coincidir con el número de filas; las incertidumbres deben ser estrictamente positivas. Las unidades de radiancia son las heredadas del caso Fortran; **la implementación no convierte automáticamente una unidad arbitraria a otra**.

## 5. `resolution.txt`

Primera línea: número de filas; filas siguientes: **número de onda** y **FWHM** instrumental, en las unidades de número de onda del modelo.

El lector ordena los números de onda, interpola linealmente la FWHM a cada punto de `.spe` y rechaza espectros que queden fuera de su cobertura. No extrapola.

La convolución instrumental también depende del `gaussian_width` del preset, de modo que debe compararse con la configuración original.

## 6. `corrch4.txt`

Archivo de coeficientes del continuo CIA H₂–CH₄, leído por `cia_python.py`. El modelo usa también H₂–H₂ y H₂–He. El lector espera el formato histórico de los coeficientes; debe mantenerse el archivo correspondiente al caso. `--check-inputs` no garantiza que el contenido de este fichero se haya validado porque el cálculo CIA ocurre durante la síntesis.

## 7. Tablas `.opa`

Por cada especie del `.fuel` se espera:

```text
<opacity-dir>/jupiter_ch4_<gas_name>.opa
```

Tablas binarias de `float32` con rejillas de 361 niveles, 5 temperaturas tabuladas por nivel y bloques espectrales de 10 cm⁻¹. El orden y la codificación son parte de la compatibilidad histórica: no se deben transformar a texto ni modificar su estructura sin adaptar `opacity.py` y la información de rejilla.

La opción `--endian` permite elegir el orden de bytes (`little` por defecto). La opción `--opacity-metadata` admite un JSON cuyas claves son nombres de gas y cuyos valores contienen campos admitidos por `OpacityGrid`: `origin`, `blocks`, `nfreq`, `dnu`, `ntemp`, `nlevel`.

Las tablas completas pueden ocupar mucho espacio y no se incluyen en el repositorio público. Se usan rutas configurables, **no rutas personales escritas en el código**.

## 8. Resultados

### `jupiter_ch4.res`

Cuatro columnas: `nu_cm^-1`, `observed`, `synthetic`, `error`.

### `jupiter_ch4.inv`

- `p_bar`, `T_K`.
- Una columna por gas recuperado con flag 2 o 3, según el orden de `.fuel`, expresado como `gas_profile_over_scale`.
- `tau_cloud` si `cloud` está activo.
- Una columna de incertidumbre `sigma_<parameter>` por parámetro invertido.

El archivo **contiene una cabecera comentada** que identifica exactamente las columnas de la ejecución; conviene usarla en scripts de análisis en vez de asumir posiciones fijas.

### `jupiter_ch4.sum`

Una sola fila con iteración final, χ² medio y `dof_T`, `dof_<gas>`, `dof_aerosol` según los parámetros realmente invertidos.

### `jupiter_ch4.avg`

Una fila por cada nivel del 1 al 361, y una columna por parámetro con el **índice 1-based en el que la columna del averaging kernel alcanza su máximo**. Es un resumen del kernel, **no una matriz de averaging kernels completa**.

### `diagnostics.json`

Incluye `chi2_mean`, `iterations`, `parameters`, `correlation_length`, `filters`, `molecular_jacobian` y `history`, con diagnóstico por iteración. El campo `relative_chi2_change` es firmado; su primera aparición es `null` (presentado como `N/A` en pantalla).

## 9. Qué comprobar antes de comparar Fortran y Python

Verificar **en ambas ejecuciones** el caso, número de puntos, intervalo espectral, especies y flags, perfil a priori, factores de escala, preset, parámetros `c_*` y `factor`, filtro, jacobianos, criterio de parada y versión del programa. Un χ² parecido **no demuestra** que las distribuciones verticales deban coincidir si las regularizaciones son distintas.
