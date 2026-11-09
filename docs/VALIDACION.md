# Estado de validación científica — Fortran frente a Python

**Código documentado:** Jupiter RT Python v0.4.3.  
**Corte:** 9 de noviembre de 2026.

## 1. Principios de comparación

La reproducibilidad se evalúa por separado en:

1. **Modelo directo:** espectro sintético calculado a partir de una atmósfera.
2. **Inversión:** perfiles finales y trayectoria iterativa.
3. **Contenido de información:** DOFS, errores estimados y averaging kernels.
4. **Configuración:** idéntico caso, espectro, instrumento, `c_*`, `factor`, filtros y criterio de parada.

Un ajuste espectral similar no prueba por sí solo que la inversa sea única. En particular, los valores de regularización distintos pueden producir soluciones verticales diferentes con χ² casi idéntico.

Los datos `.opa`, observaciones y salidas extensas de Fortran **no se distribuyen en este paquete de código**. Los números siguientes son referencias de comparaciones efectuadas durante el desarrollo; sin esos archivos no son pruebas automáticas ejecutables por un tercero.

## 2. Continuo CIA

El modelo de continuo Python implementa H₂–H₂, H₂–He y H₂–CH₄ sin ejecutar el Fortran original. Se compararon componentes frente a resultados generados por el código Fortran de referencia.

| Componente | Mayor diferencia relativa observada |
|---|---:|
| H₂–H₂ | 5,644 × 10⁻⁷ |
| H₂–He | 3,322 × 10⁻⁷ |
| H₂–CH₄ | 3,411 × 10⁻⁷ |

El ensayo cubrió temperaturas de 50 a 757,1 K y frecuencias entre 500 y 2000 cm⁻¹, en el marco de las pruebas iniciales de la traducción. Una comparación adicional de profundidad óptica de continuo en niveles atmosféricos seleccionados produjo diferencias relativas inferiores a 2,7 × 10⁻⁷.

**Alcance:** confirma concordancia numérica de estas rutinas en el conjunto probado; no certifica por sí solo las regiones de extrapolación física fuera de los intervalos de validez declarados en las aproximaciones originales.

## 3. Retrieval de temperatura

Caso histórico con **361 niveles y 373 muestras espectrales**:

| Indicador | Fortran | Python |
|---|---:|---:|
| χ² medio final | ≈ 3150,2481 | ≈ 3149,9217 |
| Diferencia RMS en temperatura | \- | ≈ 0,028 K |

La discrepancia espectral y térmica entre ambas implementaciones es pequeña para el ensayo estudiado. Se utilizó el tratamiento térmico distribuido entre niveles adyacentes, que no es una transcripción literal de todos los detalles del jacobiano original.

## 4. Abundancia molecular normal: C₂H₂ (`c2h3`)

Caso validado con temperatura troposférica y abundancia, **361 niveles y 133 puntos** entre 685,13 y 719,93 cm⁻¹:

| Indicador | Fortran | Python |
|---|---:|---:|
| χ² medio | ≈ 0,229879 | ≈ 0,229197 |
| Iteración final | 9 | 9 |
| Diferencia RMS temperatura | \- | ≈ 0,0173 K |
| Diferencia RMS relativa abundancia | \- | ≈ 0,128 % |

Esta comparación utilizó la configuración histórica para abundancia **`c_temperature=3`, `c_molecule=3`, `factor=3`**, con `--jacobian original`. La etiqueta del gas y del fichero de opacidad se conserva como `c2h3`.

**Conclusión:** excelente concordancia numérica en el caso probado, sin que ello valide de forma independiente la fórmula heredada del jacobiano molecular.

## 5. Estudio de los jacobianos moleculares

Se compararon `original`, `analytical` y `finite_difference` en el caso de abundancia normal:

| Método | χ² medio final | Iteración final | DOFS T | DOFS gas |
|---|---:|---:|---:|---:|
| `original` | ≈ 0,229197 | 9 | ≈ 1,01186 | ≈ 1,41289 |
| `analytical` | ≈ 0,231308 | 8 | ≈ 1,00608 | ≈ 1,41258 |
| `finite_difference` | ≈ 0,231308 | 8 | ≈ 1,00608 | ≈ 1,41258 |

La divergencia final de χ² refleja en parte el **diferente punto de parada**: en la iteración 8 el método original tenía χ² ≈ 0,231327, prácticamente igual al analítico.

`analytical` y `finite_difference` coincidieron con precisión muy elevada. Pero la segunda alternativa diferencia de forma centrada **una fórmula de transmisión** y comparte código con la primera; **todavía no se ha hecho una perturbación independiente del forward completo** para validar su exactitud frente a la derivada física. Por otra parte, `original` emplea una rejilla espectral submuestreada y los nuevos métodos la completa, lo que puede contribuir a las discrepancias.

La distribución vertical de máximos en `jupiter_ch4.avg` no es idéntica, especialmente para la molécula. Los DOFS cercanos no implican kernels verticales iguales; para afirmarlo sería necesario guardar y comparar las matrices completas de averaging kernels.

**Decisión actual:** mantener `original` como jacobiano molecular predeterminado.

## 6. C₆H₆: caso especial comprobado

Se comparó el caso correcto `00-01` de Fortran con la ejecución Python correspondiente:

- **361 niveles; 54 puntos** entre 671,11 y 676,94 cm⁻¹.
- Parámetros: temperatura y C₆H₆.
- Ambas ejecuciones terminaron en la **iteración 8**.

| Indicador | Fortran | Python |
|---|---:|---:|
| χ² medio final | ≈ 0,0803665 | ≈ 0,080999 |
| DOFS temperatura | ≈ 0,954 | ≈ 0,95416 |
| DOFS C₆H₆ | ≈ 0,953 | ≈ 0,95273 |

La diferencia relativa del χ² medio es aproximadamente **0,79 %**. La máxima diferencia entre espectros sintéticos fue de aproximadamente **1,32 × 10⁻¹¹**, equivalente a **0,0085 veces** la incertidumbre observacional correspondiente. La temperatura difirió RMS ≈ **0,216 K**, con máxima diferencia ≈ **0,814 K**.

**Advertencia:** los perfiles VMR de C₆H₆ mostraron discrepancias verticales apreciables (del orden de factores de dos o tres en algunas regiones) pese al excelente acuerdo espectral. Según la evaluación del caso, las ejecuciones usaron **valores de regularización diferentes**, por lo que no hay base para atribuir esas diferencias a un fallo del código. El equipo decidió **conservar la configuración 0.4.3 sin cambios**.

Esta prueba respalda el funcionamiento del modelo radiativo y de la vía de inversión especial de C₆H₆; **no debe citarse como reproducción exacta del perfil Fortran bajo parámetros matemáticos idénticos**.

### Escalado especial

La implementación lee el VMR de C₆H₆ de la **última columna** del `.pta`, aplica los factores comunes y **multiplica adicionalmente por 50 solo para `c6h6`**, reproduciendo un escalado histórico identificado en su Fortran. No se aplica a HCN, CO₂ ni H₂O.

## 7. Aerosoles y restantes moléculas especiales

| Componente | Situación al cierre | Qué queda |
|---|---|---|
| Aerosoles | Implementación heredada preservada, con perfil inicial cero | Comparación completa del retrieval con referencia Fortran y verificación de filtros y kernels. |
| CO₂ | **Siguiente validación prevista** | Correr el caso correcto y comparar χ², espectro, perfil, DOFS e incertidumbres. |
| HCN | No validado individualmente | Revisar su rejilla `.opa`, filtro específico y caso Fortran. |
| H₂O | No validado individualmente | Revisar su rejilla `.opa`, filtro específico y caso Fortran. |

El plan actual es utilizar **CO₂ como caso representativo** de la infraestructura compartida por HCN, CO₂ y H₂O. Un resultado favorable de CO₂ ayudaría a validar el flujo compartido, pero **no demostraría automáticamente** la corrección de cada opacidad específica.

## 8. Lista de comprobación para el próximo caso (CO₂)

- [ ] Asegurar que `.fuel` selecciona `T` y `co2` y que `co2` es el último gas.
- [ ] Confirmar que `.pta` contiene solo presión, T y una columna VMR por gas; no incluir aerosoles.
- [ ] Confirmar número e intervalo de puntos espectrales con `.spe`.
- [ ] Registrar el preset `co2` (`gaussian_width=1.2`, `c_special=100`, `factor=3`) y `--jacobian original`.
- [ ] Registrar regularización y configuración de Fortran; igualarlas si se pretende una comparación numérica uno a uno.
- [ ] Comparar espectros en unidades de error observacional, no solo χ².
- [ ] Comparar perfil inicial y final, incluyendo factores de escala del `.fuel`.
- [ ] Comparar DOFS y máximos de averaging kernels (sin equipararlos a matrices completas).
- [ ] Mantener las diferencias pendientes documentadas antes de declarar validación total.

## 9. Límites de las afirmaciones de validación

Esta documentación contiene comparaciones de casos específicos, no un certificado de reproducibilidad para todas las latitudes, moléculas, resoluciones instrumentales ni rangos termodinámicos. Se recomienda conservar referencias de versión de código y parámetros para cada ejecución y, si se publican resultados científicos, repetir las pruebas con los archivos y precisiones relevantes.
