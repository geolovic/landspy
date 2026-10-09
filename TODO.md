# Tareas pendientes

Registradas el 9 de octubre de 2026 tras la revisión de documentación.
Estas tareas quedan pendientes; no se ha autorizado todavía modificar su comportamiento.

- [ ] **HCurve — cálculos y momentos.** Corregir los índices intercambiados de asimetría/curtosis y sus variantes de densidad; incluir el último intervalo en la integración de HI2. Confirmar el cambio de resultados y añadir pruebas numéricas.
- [ ] **API — argumentos y coordenadas.** Hacer que Grid/DEM respeten `band`; aplicar la inversión de orden en `Channel.getXY(head=False)`; corregir la eliminación de duplicados en Flow/Network.snapPoints para conservar XY y columnas adicionales. Decidir el contrato de la columna de índice de destino y probarlo.
- [ ] **Archivos — persistencia.** Corregir el guardado y recuperación de los límites de regresiones de Channel; evitar pérdida de elevaciones negativas o precisión en Flow.save/load. Definir el formato nuevo y mantener, cuando sea posible, la lectura de archivos antiguos; cubrir ambos formatos con pruebas.

Las demás discrepancias detectadas permanecen registradas en [DOC_API_REVIEW.md](DOC_API_REVIEW.md). No modificar ejemplos ni tutoriales sin petición expresa.
