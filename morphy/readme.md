# Morfología matemática y análisis de imágenes

**Fundamentos teóricos y bases para su aplicación en datos georreferenciados de variables físicas**

---

## Resumen ejecutivo

Este reporte presenta una exposición sistemática de la morfología matemática y su papel dentro del análisis de imágenes, priorizando sus fundamentos matemáticos, geométricos y algorítmicos. Se desarrollan los operadores básicos y avanzados, sus propiedades estructurales y su relación con tareas de segmentación, extracción de objetos y caracterización geométrica. Posteriormente, se discute cómo este marco teórico se adapta al trabajo con datos georreferenciados de variables físicas, tales como temperatura superficial del mar, viento, precipitación, batimetría y otros campos raster geocientíficos. Finalmente, se propone un entorno computacional centrado en Python, acompañado de herramientas auxiliares para manejo geoespacial, teledetección y procesamiento especializado.

---

## Tabla de contenidos

- [Resumen](#resumen)
- [Introducción](#introducción)
- [Fundamentos matemáticos de la morfología matemática](#fundamentos-matemáticos-de-la-morfología-matemática)
  - [Conjuntos, imágenes y estructuras de orden](#conjuntos-imágenes-y-estructuras-de-orden)
  - [Imágenes binarias, en niveles de gris y multivariadas](#imágenes-binarias-en-niveles-de-gris-y-multivariadas)
  - [Elemento estructurante](#elemento-estructurante)
  - [Erosión y dilatación](#erosión-y-dilatación)
  - [Apertura y cierre](#apertura-y-cierre)
  - [Dualidad y propiedades fundamentales](#dualidad-y-propiedades-fundamentales)
  - [Gradiente morfológico, top-hat y black-hat](#gradiente-morfológico-top-hat-y-black-hat)
  - [Transformación hit-or-miss](#transformación-hit-or-miss)
  - [Reconstrucción morfológica](#reconstrucción-morfológica)
  - [Esqueletización, adelgazamiento y poda](#esqueletización-adelgazamiento-y-poda)
  - [Operadores conectados](#operadores-conectados)
  - [Granulometrías y análisis multiescala](#granulometrías-y-análisis-multiescala)
- [Relación entre morfología matemática y análisis de imágenes](#relación-entre-morfología-matemática-y-análisis-de-imágenes)
  - [La morfología dentro del análisis de imágenes](#la-morfología-dentro-del-análisis-de-imágenes)
  - [Preprocesamiento y realce de estructuras](#preprocesamiento-y-realce-de-estructuras)
  - [Segmentación y extracción de objetos](#segmentación-y-extracción-de-objetos)
  - [Watershed y morfología](#watershed-y-morfología)
  - [Diferencias frente a otros enfoques](#diferencias-frente-a-otros-enfoques)
- [Aspectos conceptuales del análisis de imágenes digitales](#aspectos-conceptuales-del-análisis-de-imágenes-digitales)
  - [Discretización y representación matricial](#discretización-y-representación-matricial)
  - [Vecindad y conectividad](#vecindad-y-conectividad)
  - [Resolución y escala](#resolución-y-escala)
  - [Ruido, contraste y umbralización](#ruido-contraste-y-umbralización)
  - [Limitaciones derivadas de la discretización](#limitaciones-derivadas-de-la-discretización)
- [Aplicación a datos georreferenciados de variables físicas](#aplicación-a-datos-georreferenciados-de-variables-físicas)
  - [Del espacio de píxeles al espacio físico](#del-espacio-de-píxeles-al-espacio-físico)
  - [Resolución espacial y anisotropía](#resolución-espacial-y-anisotropía)
  - [Máscaras, bordes y valores faltantes](#máscaras-bordes-y-valores-faltantes)
  - [Proyecciones y sistemas de referencia](#proyecciones-y-sistemas-de-referencia)
  - [Raster geofísicos frente a imágenes convencionales](#raster-geofísicos-frente-a-imágenes-convencionales)
  - [Series espacio-temporales](#series-espacio-temporales)
  - [Ejemplos de interés geocientífico](#ejemplos-de-interés-geocientífico)
- [Lineamientos metodológicos para aplicaciones futuras](#lineamientos-metodológicos-para-aplicaciones-futuras)
  - [Selección del elemento estructurante](#selección-del-elemento-estructurante)
  - [Elección de conectividad](#elección-de-conectividad)
  - [Sensibilidad a umbral y resolución](#sensibilidad-a-umbral-y-resolución)
  - [Interpretación física](#interpretación-física)
  - [Reproducibilidad](#reproducibilidad)
- [Software computacional recomendado](#software-computacional-recomendado)
  - [Python como entorno principal](#python-como-entorno-principal)
  - [Herramientas auxiliares](#herramientas-auxiliares)
  - [Comparación sintética](#comparación-sintética)
- [Conclusiones y recomendaciones para investigación aplicada](#conclusiones-y-recomendaciones-para-investigación-aplicada)
- [Apéndice: Esquema sugerido de flujo de trabajo en Python](#apéndice-esquema-sugerido-de-flujo-de-trabajo-en-python)
- [Bibliografía](#bibliografía)

---

## Resumen

La morfología matemática constituye un marco formal para el análisis de estructuras geométricas en imágenes. Su valor radica en que permite estudiar objetos, bordes, huecos, conectividad y escalas espaciales a partir de operadores no lineales definidos sobre conjuntos, funciones y retículas. A diferencia de métodos puramente espectrales o estadísticos, la morfología matemática incorpora explícitamente la noción de forma, vecindad y conectividad. Por ello, resulta especialmente útil cuando el interés científico está puesto en identificar entidades espaciales coherentes, delimitar regiones conectadas, suprimir artefactos o caracterizar patrones con interpretación física.

En este documento se prioriza primero la teoría: se introducen las bases matemáticas, las imágenes binarias y en niveles de gris, los elementos estructurantes y los principales operadores morfológicos. Después se explica la relación de estas ideas con el análisis de imágenes digitales más amplio, incluyendo segmentación, watershed, etiquetado de componentes y análisis multiescala. Solo en una segunda parte se discute la aplicación a datos georreferenciados, donde intervienen factores adicionales como resolución espacial, proyecciones, bordes, valores faltantes, máscaras y escalas físicas. Se concluye con recomendaciones metodológicas y una revisión del software computacional, destacando a Python como entorno principal.

---

## Introducción

La morfología matemática es una teoría del análisis de formas desarrollada inicialmente en el contexto del tratamiento de imágenes binarias y posteriormente extendida a imágenes en niveles de gris, multibanda y estructuras de mayor dimensionalidad. Su núcleo conceptual consiste en estudiar cómo un objeto o una imagen cambia al ser *sondeado* por una forma elemental llamada *elemento estructurante*. Esta idea convierte la geometría local y la conectividad en piezas centrales del análisis.

Desde una perspectiva amplia, el análisis de imágenes busca extraer información significativa a partir de datos espaciales discretizados. Sin embargo, no todas las preguntas sobre una imagen son equivalentes. Algunas se refieren a intensidad media, varianza o texturas; otras, en cambio, se centran en la presencia de bordes, la delimitación de regiones, el número de objetos conectados, el tamaño de huecos o la persistencia de estructuras a distintas escalas. La morfología matemática es especialmente adecuada para este segundo grupo de problemas.

En ciencias e ingeniería, este marco resulta relevante porque muchos fenómenos se manifiestan como patrones espaciales: filamentos, parches, manchas, frentes, núcleos intensos, bandas, regiones conectadas o estructuras ramificadas. Tales entidades aparecen tanto en imágenes tradicionales como en datos científicos georreferenciados. En oceanografía, meteorología, hidrología, geología y teledetección, es frecuente trabajar con campos raster de temperatura, clorofila, precipitación, topografía, reflectancia o intensidad de eventos extremos. En todos esos casos, no basta con analizar solamente valores puntuales; también es necesario comprender la organización espacial de las estructuras.

El presente reporte adopta deliberadamente un orden de exposición donde la teoría precede a la aplicación. Esta decisión es metodológicamente útil por dos razones. Primero, porque evita reducir la morfología matemática a un catálogo de funciones computacionales. Segundo, porque permite entender cómo y por qué ciertas operaciones producen determinados efectos geométricos antes de trasladarlas a contextos más complejos, como los datos georreferenciados de variables físicas.

---

## Fundamentos matemáticos de la morfología matemática

### Conjuntos, imágenes y estructuras de orden

La formulación clásica de la morfología matemática parte de la teoría de conjuntos. Si una imagen binaria se interpreta como el conjunto de píxeles pertenecientes al objeto, entonces muchas operaciones sobre la imagen pueden formularse como operaciones entre conjuntos. Sea $E$ un espacio discreto, usualmente una grilla bidimensional

```math
E \subseteq \mathbb{Z}^2.
```


Una imagen binaria puede representarse como un subconjunto $A \subseteq E$, donde cada punto de $A$ corresponde a un píxel activo u objeto.

Con el desarrollo de la teoría, este enfoque se amplió hacia imágenes en escala de grises, modeladas como funciones

```math
f : E \to \mathbb{R}.
```


En este caso, la morfología deja de actuar solo sobre pertenencia o no pertenencia, y pasa a operar sobre órdenes parciales entre funciones. Esto conduce a una formulación sobre *retículas completas*, donde existen operaciones de supremo e ínfimo. La importancia de este marco radica en que muchos operadores morfológicos pueden definirse abstractamente como aplicaciones crecientes que preservan ciertas propiedades de orden.

#### Retículas y operadores

Una retícula es un conjunto parcialmente ordenado donde cada par de elementos admite un supremo y un ínfimo. En el contexto morfológico, esto permite trabajar con familias de imágenes y definir operadores en términos de máximos y mínimos locales. Conceptualmente, erosiones y dilataciones son transformaciones ligadas a ínfimos y supremos respecto de vecindades determinadas por el elemento estructurante.

### Imágenes binarias, en niveles de gris y multivariadas

#### Imágenes binarias
Una imagen binaria puede verse como una función indicadora

```math
\chi_A(x)=
\begin{cases}
1, & x \in A,\\
0, & x \notin A.
\end{cases}
```


Esta representación es útil cuando el interés está puesto en presencia/ausencia de un objeto, o cuando una imagen continua ha sido umbralizada para obtener una máscara.

#### Imágenes en niveles de gris
Una imagen en niveles de gris se representa por una función escalar $f(x)$ definida sobre la grilla. Aquí la intensidad del píxel importa, y los operadores morfológicos actúan desplazando máximos y mínimos locales bajo la acción del elemento estructurante.

#### Imágenes multivariadas
En aplicaciones reales es común trabajar con imágenes multicanal o campos multivariados, por ejemplo bandas espectrales de satélite o variables físicas asociadas. La extensión morfológica a este caso requiere definir un orden adecuado entre vectores o recurrir a estrategias por componente, transformaciones previas o enfoques vectoriales especializados.

### Elemento estructurante

El elemento estructurante es una pieza central de la teoría. Se trata de un conjunto pequeño $B \subseteq E$ o, más generalmente, de una función auxiliar que define la geometría local con la cual la imagen será interrogada. En términos intuitivos, el elemento estructurante fija la forma, escala y orientación de la exploración local.

Por ejemplo, si $B$ es un disco pequeño, el operador tenderá a preservar o enfatizar estructuras aproximadamente isotrópicas. Si $B$ es una línea orientada, el análisis será sensible a estructuras alargadas en una dirección concreta. Esta relación entre la forma del elemento estructurante y la interpretación geométrica del resultado es una de las razones por las cuales la morfología matemática es especialmente valiosa en problemas con significado físico.

### Erosión y dilatación

#### Erosión en imágenes binarias
Sea $A \subseteq E$ una imagen binaria y $B$ un elemento estructurante. La erosión de $A$ por $B$ se define como

```math
A \ominus B = \{x \in E : B_x \subseteq A\},
```


donde $B_x$ denota la traslación de $B$ centrada en $x$.

Interpretación: un punto permanece en la imagen erosionada solo si el elemento estructurante cabe completamente dentro del objeto al colocarse en dicho punto. La erosión tiende a contraer los objetos, eliminar protrusiones pequeñas y romper conexiones delgadas.

#### Dilatación en imágenes binarias
La dilatación de $A$ por $B$ se define como

```math
A \oplus B = \{x \in E : (\hat{B})_x \cap A \neq \varnothing\},
```


donde $\hat{B}$ es la reflexión de $B$ respecto del origen.

Interpretación: la dilatación expande el objeto, rellena huecos pequeños y conecta regiones próximas si la distancia entre ellas es menor que la escala impuesta por $B$.

#### Versión en niveles de gris
Para una imagen en niveles de gris $f$ y un elemento estructurante plano $B$, la dilatación y la erosión pueden escribirse como

```math
(\delta_B f)(x) = \sup_{b \in B} f(x-b),
\qquad
(\varepsilon_B f)(x) = \inf_{b \in B} f(x+b).
```


Estas expresiones muestran que la dilatación actúa como un máximo local móvil y la erosión como un mínimo local móvil.

### Apertura y cierre

La apertura se define como una erosión seguida de una dilatación:

```math
A \circ B = (A \ominus B)\oplus B.
```


La apertura es *anti-extensiva*, es decir,

```math
A \circ B \subseteq A.
```


Geométricamente, elimina detalles brillantes o salientes pequeños que no pueden contener al elemento estructurante. Además, suaviza contornos y rompe puentes estrechos.

El cierre se define como una dilatación seguida de una erosión:

```math
A \bullet B = (A \oplus B)\ominus B.
```


El cierre es *extensivo*,

```math
A \subseteq A \bullet B.
```


Tiende a rellenar huecos oscuros pequeños, cerrar grietas y conectar discontinuidades estrechas.

### Dualidad y propiedades fundamentales

Una parte importante de la elegancia de la teoría morfológica es la dualidad entre operadores. Bajo complemento y reflexión del elemento estructurante, erosión y dilatación son duales; lo mismo ocurre con apertura y cierre. Esta dualidad no es solo formal: expresa que los operadores pueden actuar ya sea sobre el objeto o sobre el fondo, según la interpretación del problema.

Entre las propiedades más relevantes destacan:

- **Crecimiento:** si $A \subseteq C$, entonces
  
```math
A \ominus B \subseteq C \ominus B, \qquad A \oplus B \subseteq C \oplus B.
```

- **Anti-extensividad de la apertura:**
  
```math
A \circ B \subseteq A.
```

- **Extensividad del cierre:**
  
```math
A \subseteq A \bullet B.
```

- **Idempotencia:**
  
```math
(A \circ B)\circ B = A \circ B, \qquad (A \bullet B)\bullet B = A \bullet B.
```


Estas propiedades son esenciales para interpretar los operadores como filtros geométricos estables.

### Gradiente morfológico, top-hat y black-hat

El gradiente morfológico es una medida de contraste espacial construida a partir de dilatación y erosión:

```math
\nabla_B f = \delta_B f - \varepsilon_B f.
```


Este operador resalta transiciones locales y, por ello, suele emplearse en detección de bordes.

El *white top-hat* se define como

```math
\mathrm{WTH}_B(f)=f-(f\circ B),
```


y destaca detalles brillantes pequeños respecto de la escala impuesta por $B$.

El *black-hat* se define como

```math
\mathrm{BTH}_B(f)=(f\bullet B)-f,
```


y enfatiza detalles oscuros pequeños.

Ambos operadores son muy útiles cuando se desea separar estructuras locales de un fondo suavizado.

### Transformación hit-or-miss

La transformación *hit-or-miss* está diseñada para detectar configuraciones locales específicas en imágenes binarias. Emplea dos elementos estructurantes: uno para el objeto y otro para el fondo. De manera esquemática, busca coincidencias simultáneas entre una forma objetivo y su entorno. Esta idea es importante en reconocimiento de patrones simples, detección de esquinas, terminaciones o configuraciones topológicas particulares.

### Reconstrucción morfológica

La reconstrucción morfológica es un procedimiento más selectivo que las operaciones básicas. Parte de una *imagen marcador* y una *imagen máscara*. Mediante dilataciones geodésicas iteradas restringidas por la máscara, reconstruye componentes conectadas compatibles con el marcador. Conceptualmente, esto permite eliminar artefactos preservando mejor la forma global de las estructuras relevantes.

Su importancia práctica es grande porque muchas aperturas y cierres clásicos deforman geometrías de manera no deseada. En contraste, la reconstrucción tiende a ser más fiel a la conectividad y estructura de los objetos que sobreviven al criterio impuesto.

### Esqueletización, adelgazamiento y poda

La esqueletización busca representar un objeto mediante un conjunto más delgado que preserve, idealmente, su topología esencial. El esqueleto intenta condensar la forma en una estructura central, útil para describir ramificaciones, ejes o conectividad.

El adelgazamiento es un procedimiento iterativo que remueve píxeles del borde sin alterar la conectividad global. La poda, por su parte, elimina ramas espurias o demasiado cortas del resultado esqueletizado, reduciendo ruido geométrico.

Estos operadores son especialmente útiles cuando el interés se centra en redes, trayectorias, filamentos o geometrías lineales.

### Operadores conectados

Los operadores conectados forman una clase avanzada de filtros morfológicos que actúan sobre componentes conectadas más que sobre ventanas locales fijas. Su ventaja principal es que permiten eliminar estructuras según atributos como área, volumen, altura o persistencia, preservando contornos de componentes que no son eliminadas.

Esto los vuelve conceptualmente muy atractivos para aplicaciones científicas donde interesa suprimir objetos irrelevantes por tamaño o intensidad sin distorsionar la geometría de los restantes.

### Granulometrías y análisis multiescala

Una granulometría es una familia de aperturas parametrizadas por el tamaño del elemento estructurante:

```math
\{\gamma_{\lambda}\}_{\lambda \ge 0},
```


donde $\lambda$ representa la escala. El análisis de cómo la imagen cambia con $\lambda$ permite inferir distribuciones de tamaño de estructuras. Esta idea es central en análisis multiescala: en vez de estudiar la imagen a una sola resolución geométrica, se examina qué rasgos persisten, emergen o desaparecen a distintas escalas.

En problemas geofísicos, esta perspectiva es particularmente útil porque los fenómenos pueden exhibir una jerarquía espacial clara: parches pequeños embebidos en regiones grandes, frentes delgados sobre campos amplios o estructuras coherentes cuya interpretación depende de la escala física seleccionada.

---

## Relación entre morfología matemática y análisis de imágenes

### La morfología dentro del análisis de imágenes

El análisis de imágenes comprende adquisición, preprocesamiento, realce, segmentación, extracción de características, clasificación e interpretación. La morfología matemática se inserta en varias de estas etapas, pero destaca especialmente en aquellas relacionadas con estructura geométrica y conectividad.

A diferencia de técnicas puramente lineales, como ciertas convoluciones o filtros de Fourier, los operadores morfológicos son no lineales y están guiados por la forma. Esto les permite capturar propiedades topológicas y geométricas que serían menos transparentes en un enfoque exclusivamente espectral.

### Preprocesamiento y realce de estructuras

Antes de segmentar o medir una imagen, suele ser necesario suprimir ruido, rellenar discontinuidades o destacar objetos pequeños. La apertura y el cierre cumplen aquí un papel central: uno elimina protuberancias y ruido brillante pequeño; el otro rellena huecos oscuros y cierra fracturas estrechas. El gradiente morfológico, por su parte, es útil para resaltar contornos.

### Segmentación y extracción de objetos

La segmentación consiste en separar regiones significativas del fondo o de otras regiones. La morfología matemática participa en este proceso de varias maneras:

- refinando máscaras obtenidas por umbralización;
- separando objetos que se tocan;
- conectando regiones fragmentadas;
- eliminando componentes pequeñas irrelevantes;
- generando marcadores para segmentación basada en watershed.

Una vez segmentadas las regiones, el etiquetado de componentes conectados permite identificar objetos individuales y calcular atributos como área, perímetro, elongación, compacidad u orientación.

### Watershed y morfología

El algoritmo *watershed* interpreta una imagen como una superficie topográfica. Los mínimos locales actúan como cuencas, y las fronteras entre cuencas delimitan regiones. La conexión con la morfología es profunda: el watershed puede formularse morfológicamente y combinarse con reconstrucción, marcadores y gradientes morfológicos. En la práctica, esto resulta muy útil para separar objetos cercanos o superpuestos.

### Diferencias frente a otros enfoques

La morfología matemática no reemplaza otros marcos de análisis, sino que los complementa. Frente a enfoques estadísticos, aporta una representación explícita de forma y conectividad. Frente a métodos espectrales, enfatiza estructura local y geometría. Frente a modelos de aprendizaje automático, ofrece interpretabilidad directa y control explícito sobre la escala y la forma del análisis. Esta combinación la vuelve muy valiosa como capa metodológica intermedia entre datos y significado físico.

---

## Aspectos conceptuales del análisis de imágenes digitales

### Discretización y representación matricial

Toda imagen digital es una aproximación discreta de un campo continuo. Si el fenómeno físico subyacente vive en un dominio continuo, la imagen observada solo registra valores sobre una grilla finita. Esto implica pérdida de información, aliasing potencial y dependencia respecto de la resolución espacial.

Matemáticamente, una imagen puede representarse por una matriz

```math
I = (I_{ij})_{i=1,\dots,m;\,j=1,\dots,n},
```


donde cada entrada corresponde a un píxel. En contextos multibanda o espacio-temporales, esta representación se extiende a tensores.

### Vecindad y conectividad

La conectividad define qué significa que dos píxeles estén "conectados". En 2D, las nociones más comunes son conectividad-4 y conectividad-8. En 3D aparecen variantes como conectividad-6, 18 o 26. Esta decisión no es menor: afecta el conteo de componentes, la presencia de puentes entre regiones y la topología inferida a partir de la imagen.

### Resolución y escala

La resolución espacial fija el tamaño del elemento más pequeño potencialmente detectable. Una estructura de escala subpíxel no podrá reconstruirse fielmente. Además, la elección del elemento estructurante debe interpretarse en relación con la resolución: un disco de radio 3 píxeles tiene significados físicos distintos según el tamaño real de píxel.

### Ruido, contraste y umbralización

El ruido afecta la estabilidad de la segmentación y del análisis morfológico. Pequeñas fluctuaciones de intensidad pueden generar objetos espurios o conexiones artificiales. Por ello, la combinación entre umbralización, filtrado previo y morfología es una práctica habitual. Sin embargo, debe hacerse con cuidado: un umbral inadecuado puede crear artefactos topológicos que luego la morfología solo maquilla, en lugar de resolver de manera conceptualmente correcta.

### Limitaciones derivadas de la discretización

La discretización puede distorsionar ángulos, curvaturas, orientaciones y tamaños. También puede introducir anisotropía cuando los píxeles no son cuadrados o cuando las dimensiones tienen escalas distintas. Estas limitaciones son especialmente relevantes al transferir operadores concebidos geométricamente a grillas reales de observación.

---

## Aplicación a datos georreferenciados de variables físicas

### Del espacio de píxeles al espacio físico

En una imagen convencional, el píxel suele ser una unidad abstracta de muestreo. En cambio, en un raster georreferenciado cada celda representa una porción concreta del espacio físico. Por ello, el análisis morfológico no puede quedarse en la escala de píxel: debe reinterpretarse en términos de kilómetros, metros o grados, según el sistema de referencia.

Esto implica que la elección del elemento estructurante no debe hacerse solo por conveniencia computacional. Debe responder a una escala física del fenómeno bajo estudio. Por ejemplo, detectar parches térmicos costeros, filamentos de clorofila o núcleos de precipitación requiere elementos estructurantes coherentes con el tamaño esperado de las estructuras.

### Resolución espacial y anisotropía

En datos georreferenciados, la resolución puede variar con la latitud, con la proyección o entre ejes coordenados. Una celda puede no representar la misma distancia en las direcciones zonal y meridional, y ello afecta directamente la interpretación geométrica de vecindad y distancia. En estos casos, un elemento estructurante isotrópico en coordenadas de píxel puede no ser isotrópico en el espacio físico.

### Máscaras, bordes y valores faltantes

Los campos geofísicos contienen con frecuencia valores faltantes, regiones enmascaradas o bordes irregulares. En oceanografía, por ejemplo, la máscara tierra-mar introduce discontinuidades espaciales no triviales. Aplicar operadores morfológicos sin considerar estas restricciones puede producir conexiones artificiales entre zonas separadas por tierra o interpretar como fondo lo que simplemente es ausencia de dato.

Por ello, una implementación rigurosa debe decidir explícitamente cómo tratar:
- valores `NaN` o `nodata`;
- bordes del dominio;
- regiones enmascaradas;
- discontinuidades físicas o geométricas.

### Proyecciones y sistemas de referencia

La georreferenciación introduce otro nivel de complejidad: las distancias, áreas y orientaciones dependen del sistema de coordenadas y de la proyección cartográfica. Si la operación morfológica debe tener significado geométrico en unidades físicas, es necesario verificar si el dominio de trabajo preserva localmente esas magnitudes de manera adecuada. En muchos casos conviene reproyectar a una malla métrica antes de aplicar operadores sensibles a escala y forma.

### Raster geofísicos frente a imágenes convencionales

Aunque un campo de temperatura superficial del mar pueda representarse como una imagen, su interpretación difiere de la de una fotografía. En imágenes físicas, los niveles de gris no representan necesariamente reflectancia visual, sino magnitudes continuas con unidades, errores de medición, ruido instrumental, escalas de correlación y restricciones dinámicas. Así, una apertura o un cierre no solo modifica una textura visual; altera potencialmente la delimitación de un fenómeno físico.

### Series espacio-temporales

Muchas aplicaciones geocientíficas involucran secuencias temporales de mapas. En ese caso, la noción de conectividad puede extenderse al dominio espacio-temporal. Esto permite identificar eventos persistentes, objetos móviles o estructuras conectadas a través del tiempo. Sin embargo, dicha extensión exige cuidado: la escala temporal debe elegirse con el mismo rigor que la espacial, y la conectividad entre tiempos sucesivos debe tener interpretación física.

### Ejemplos de interés geocientífico

Sin entrar todavía en un estudio de caso detallado, la morfología matemática puede apoyar tareas como:
- delimitación de frentes y bordes intensos;
- detección de parches anómalos;
- limpieza de máscaras de eventos extremos;
- identificación de regiones conectadas en campos umbralizados;
- extracción de filamentos, lenguas o estructuras alargadas;
- análisis multiescala de heterogeneidad espacial;
- posprocesamiento de resultados de segmentación o clasificación.

---

## Lineamientos metodológicos para aplicaciones futuras

### Selección del elemento estructurante

La forma y tamaño del elemento estructurante deben responder al problema científico. Un disco es razonable cuando se espera isotropía local; una línea orientada puede ser más adecuada para filamentos o frentes; una cruz o una vecindad Manhattan puede tener sentido cuando la conectividad de la grilla es prioritaria. En todos los casos, el tamaño debe justificarse en unidades físicas.

### Elección de conectividad

La conectividad no es un detalle técnico menor. Puede cambiar el número de objetos detectados, la identificación de puentes y la topología de la segmentación. Por ello, la conectividad elegida debe declararse explícitamente y, de ser posible, compararse con alternativas para evaluar sensibilidad.

### Sensibilidad a umbral y resolución

Muchos análisis morfológicos se aplican sobre máscaras binarias obtenidas mediante umbralización. Por tanto, la calidad del resultado depende tanto del operador morfológico como del criterio de binarización. Es recomendable explorar sensibilidad respecto del umbral y la resolución espacial, idealmente mediante experimentos controlados.

### Interpretación física

No todo objeto morfológicamente detectado corresponde a una estructura físicamente significativa. Algunas regiones pueden surgir por ruido, remuestreo, huecos de datos o elecciones arbitrarias de preprocesamiento. Por ello, el análisis morfológico debe integrarse con conocimiento del sistema físico, no sustituirlo.

### Reproducibilidad

Un flujo de trabajo reproducible debería documentar:
- fuente y formato de los datos;
- proyección y resolución;
- tratamiento de máscaras y `nodata`;
- umbral(es) empleados;
- conectividad utilizada;
- forma y tamaño del elemento estructurante;
- software y librerías;
- criterios de validación.

---

## Software computacional recomendado

### Python como entorno principal

Python es actualmente un entorno especialmente adecuado para implementar flujos de trabajo reproducibles de morfología matemática y análisis de imágenes sobre datos científicos y georreferenciados. Su fortaleza no reside en una sola biblioteca, sino en el ecosistema articulado entre procesamiento numérico, manejo de metadatos, análisis raster y visualización.

#### NumPy
**Rol principal:** estructura base para arreglos multidimensionales y operaciones vectorizadas.

NumPy proporciona la representación matricial y tensorial sobre la cual descansan muchas operaciones. Es la base para manipular imágenes como arreglos, construir máscaras, aplicar operaciones booleanas y ejecutar transformaciones elementales de manera eficiente.

#### SciPy
**Rol principal:** rutinas científicas generales, filtrado, interpolación y operaciones numéricas complementarias.

SciPy amplía a NumPy con herramientas útiles para procesamiento de señales e imágenes, estructuras espaciales y análisis numérico. Aunque la morfología puede implementarse con bibliotecas más especializadas, SciPy sigue siendo un soporte importante dentro del flujo general.

#### scikit-image
**Rol principal:** procesamiento de imágenes científicas con énfasis en segmentación, morfología y extracción de características.

Es una de las bibliotecas más adecuadas para estudiar e implementar morfología matemática en Python. Incluye erosión, dilatación, apertura, cierre, reconstrucción, esqueletización, watershed, etiquetado y múltiples herramientas de medición geométrica.

#### OpenCV
**Rol principal:** visión por computador y operaciones de imágenes de alto rendimiento.

OpenCV ofrece implementaciones robustas de operadores morfológicos y es útil cuando se requieren pipelines rápidos, integración con otras tareas de visión o procesamiento de grandes lotes de imágenes. En aplicaciones geocientíficas debe usarse con cuidado respecto a georreferenciación y metadatos, ya que su enfoque central no es geoespacial.

#### xarray
**Rol principal:** manejo de arreglos etiquetados multidimensionales, especialmente útil para datos científicos con coordenadas y tiempo.

Para variables físicas georreferenciadas, xarray es especialmente valioso porque permite trabajar con dimensiones nombradas, coordenadas explícitas y estructuras tipo `DataArray` o `Dataset`. Esto facilita operar sobre campos espacio-temporales preservando semántica científica.

#### rioxarray
**Rol principal:** extensión geoespacial de xarray para raster con información de referencia espacial.

rioxarray permite trabajar con coordenadas de referencia, sistemas CRS y operaciones geoespaciales dentro de un flujo basado en xarray, lo que resulta muy útil cuando se combinan datos científicos y georreferenciación.

#### rasterio
**Rol principal:** lectura, escritura y manejo de raster geoespaciales.

rasterio es especialmente útil para abrir GeoTIFF, inspeccionar metadatos espaciales, manejar transformaciones afines, recortes, remuestreo y escritura de productos procesados.

#### geopandas
**Rol principal:** manipulación de datos vectoriales geoespaciales.

Aunque la morfología aquí se centra en raster, geopandas es útil para integrar resultados con líneas de costa, polígonos administrativos, regiones de estudio o capas auxiliares.

#### matplotlib
**Rol principal:** visualización científica reproducible.

La correcta visualización de máscaras, regiones conectadas, bordes y resultados multiescala es parte integral del análisis. matplotlib proporciona una base flexible para ello.

#### dask
**Rol principal:** cómputo paralelo y manejo de arreglos más grandes que memoria.

Cuando los datos son muy voluminosos, como series espacio-temporales extensas o mosaicos satelitales, dask permite escalar el procesamiento sin abandonar el ecosistema Python.

### Herramientas auxiliares

#### GDAL
GDAL es una herramienta fundamental para traducción, conversión, reproyección, submuestreo y manejo general de formatos geoespaciales raster y vectoriales. Aunque no es el eje del análisis morfológico teórico, sí es un apoyo muy importante en la preparación y transformación de datos.

#### QGIS
QGIS es útil para inspección visual, validación espacial, exploración de capas, construcción de productos intermedios y verificación cualitativa de resultados. Puede complementar bien un flujo reproducible en Python, pero no conviene que reemplace el procesamiento principal si el objetivo es automatización y trazabilidad científica.

#### Orfeo Toolbox
Orfeo Toolbox resulta especialmente valioso en teledetección. Ofrece herramientas especializadas para clasificación y perfiles morfológicos, por lo que puede ser muy útil cuando el trabajo se acerca a procesamiento avanzado de imágenes satelitales de alta resolución.

#### ITK
ITK proporciona un marco más amplio para imágenes científicas, segmentación y morfología avanzada. Puede ser conveniente en aplicaciones exigentes o cuando se requieren métodos especializados más allá del uso rutinario en Python.

### Comparación sintética

| Herramienta       | Rol principal                                               | Uso recomendado                                                                 |
|-------------------|-------------------------------------------------------------|---------------------------------------------------------------------------------|
| NumPy             | Arreglos y operaciones vectorizadas                         | Base numérica general del flujo de trabajo.                                     |
| SciPy             | Rutinas científicas y soporte numérico                      | Filtros, interpolación, operaciones auxiliares.                                 |
| scikit-image      | Morfología, segmentación y medición                         | Biblioteca principal para morfología matemática en Python.                      |
| OpenCV            | Procesamiento rápido de imágenes                            | Operaciones morfológicas eficientes, con menor énfasis geoespacial.             |
| xarray            | Arreglos etiquetados multidimensionales                     | Datos científicos con coordenadas, tiempo y metadatos.                          |
| rioxarray         | Extensión geoespacial de xarray                             | Raster georreferenciados con CRS y operaciones espaciales.                      |
| rasterio          | Entrada/salida raster                                       | GeoTIFF, metadatos espaciales, remuestreo y escritura.                          |
| geopandas         | Datos vectoriales geográficos                               | Integración con capas vectoriales auxiliares.                                   |
| matplotlib        | Visualización científica                                    | Figuras reproducibles y análisis exploratorio.                                  |
| dask              | Escalamiento computacional                                  | Datos grandes o procesamiento distribuido.                                      |
| GDAL              | Infraestructura geoespacial                                 | Conversión, reproyección y preparación de datos.                                |
| QGIS              | Exploración visual y validación                             | Inspección de resultados y control espacial.                                    |
| Orfeo Toolbox     | Teledetección especializada                                 | Perfiles morfológicos y procesamiento avanzado de imágenes.                     |
| ITK               | Imágenes científicas avanzadas                              | Segmentación y morfología especializada.                                        |

---

## Conclusiones y recomendaciones para investigación aplicada

La morfología matemática constituye un marco conceptual sólido para estudiar forma, conectividad, borde, tamaño y estructura en imágenes. Su fuerza reside en que convierte preguntas geométricas y topológicas en operadores bien definidos, con propiedades formales claras y una interpretación intuitiva potente. Por ello, no debe entenderse solo como un conjunto de técnicas prácticas, sino como una teoría del análisis estructural de imágenes.

En este reporte se ha priorizado primero esa base teórica: conjuntos, retículas, erosión, dilatación, apertura, cierre, reconstrucción, esqueletización, operadores conectados y análisis multiescala. Posteriormente, se mostró que al pasar a datos georreferenciados de variables físicas aparecen exigencias adicionales: escala física, proyecciones, anisotropía, máscaras, bordes, valores faltantes y significado físico de los objetos detectados.

Desde el punto de vista metodológico, una aplicación rigurosa exige justificar explícitamente la elección del elemento estructurante, la conectividad, la resolución y el criterio de segmentación. La reproducibilidad requiere además documentar software, parámetros y tratamiento de metadatos espaciales.

En términos computacionales, Python ofrece actualmente un ecosistema muy favorable para articular la teoría con la práctica. scikit-image proporciona una base excelente para morfología y segmentación; xarray y rioxarray facilitan el trabajo con datos científicos y georreferenciados; rasterio y GDAL permiten una interacción robusta con formatos espaciales; y QGIS, Orfeo Toolbox e ITK funcionan como herramientas auxiliares valiosas en flujos más amplios.

Como recomendación general para investigación aplicada, conviene que la morfología matemática no se use como posprocesamiento automático descontextualizado, sino como una herramienta de modelado geométrico coherente con la física del problema. Allí radica su mayor potencial en geociencias, oceanografía, meteorología, teledetección y análisis espacial avanzado.

---

## Apéndice: Esquema sugerido de flujo de trabajo en Python

A continuación se presenta un esquema conceptual de flujo de trabajo para futuras aplicaciones:

1. Lectura del dato georreferenciado (`xarray`, `rioxarray`, `rasterio`).
2. Inspección de metadatos espaciales y temporales.
3. Tratamiento de máscara, `nodata` y bordes.
4. Reproyección o remuestreo si la escala física del análisis lo exige.
5. Definición del criterio de segmentación o umbralización.
6. Elección del elemento estructurante en unidades físicamente interpretables.
7. Aplicación de operadores morfológicos.
8. Etiquetado de componentes conectados y extracción de atributos.
9. Validación visual y cuantitativa.
10. Exportación reproducible de resultados.

Ejemplo ilustrativo de pseudocódigo:

```python
import xarray as xr
import rioxarray
import numpy as np
from skimage.morphology import disk, opening, closing
from skimage.measure import label, regionprops

da = xr.open_dataset("dato.nc")["variable"]
mask = np.isfinite(da.values)
img = np.where(mask, da.values, np.nan)

# Ejemplo conceptual: umbralización
binary = np.where(img > umbral, 1, 0).astype(np.uint8)

# Elemento estructurante
se = disk(3)

# Operaciones morfológicas
binary2 = closing(opening(binary, se), se)

# Componentes conectados
lab = label(binary2, connectivity=2)
props = regionprops(lab)
```

---

## Bibliografía

1. Matheron, G. (1975). *Random Sets and Integral Geometry*. Wiley.
2. Serra, J. (1982). *Image Analysis and Mathematical Morphology*. Academic Press.
3. Serra, J. (1988). *Image Analysis and Mathematical Morphology, Volume 2: Theoretical Advances*. Academic Press.
4. Soille, P. (2003). *Morphological Image Analysis: Principles and Applications*. 2nd ed., Springer.
5. Heijmans, H. J. A. M. (1994). *Morphological Image Operators*. Academic Press.
6. Dougherty, E. R. (Ed.). (1992). *An Introduction to Morphological Image Processing*. SPIE Optical Engineering Press.
7. Vincent, L. (1993). Morphological grayscale reconstruction in image analysis: applications and efficient algorithms. *IEEE Transactions on Image Processing*, 2(2), 176–201.
8. Beucher, S., & Meyer, F. (1992). The morphological approach to segmentation: the watershed transformation. En *Mathematical Morphology in Image Processing*, 433–481. Marcel Dekker.
9. Najman, L., & Talbot, H. (Eds.). (2013). *Mathematical Morphology: From Theory to Applications*. ISTE / Wiley.
10. Najman, L., & Couprie, M. (2006). Building the component tree in quasi-linear time. *IEEE Transactions on Image Processing*, 15(11), 3531–3539.
11. Meyer, F. (1994). Topographic distance and watershed lines. *Signal Processing*, 38(1), 113–125.
12. Pesaresi, M., & Benediktsson, J. A. (2001). A new approach for the morphological segmentation of high-resolution satellite imagery. *IEEE Transactions on Geoscience and Remote Sensing*, 39(2), 309–320.
13. scikit-image developers. *scikit-image documentation*. Disponible en: https://scikit-image.org/docs/stable/
14. OpenCV developers. *OpenCV documentation*. Disponible en: https://docs.opencv.org/
15. xarray developers. *xarray documentation*. Disponible en: https://docs.xarray.dev/
16. rioxarray developers. *rioxarray documentation*. Disponible en: https://corteva.github.io/rioxarray/
17. Rasterio developers. *Rasterio documentation*. Disponible en: https://rasterio.readthedocs.io/
18. GDAL/OGR contributors. *GDAL documentation*. Disponible en: https://gdal.org/
19. QGIS Development Team. *QGIS documentation*. Disponible en: https://docs.qgis.org/
20. Orfeo Toolbox Team. *Orfeo Toolbox Cookbook*. Disponible en: https://www.orfeo-toolbox.org/CookBook/
21. Insight Software Consortium. *ITK documentation*. Disponible en: https://docs.itk.org/
