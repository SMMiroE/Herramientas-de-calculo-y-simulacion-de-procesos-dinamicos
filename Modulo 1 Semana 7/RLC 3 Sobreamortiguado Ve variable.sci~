// Curso HERRAMIENTAS DE CÁLCULO Y SIMULACIÓN DE PROCESOS DINÁMICOS
// Silvia Miró
// Modulo 1 - Semana 7
// Actividad 3: RLC Sobreamortiguado Ve variable

// --- 1. Parámetros del circuito ---
R = 2.2e3;     // Resistencia en Ohms (2.2 kOhms)
L = 10e-6;     // Inductancia en Henrys (10 microHenrys)
C = 100e-9;    // Capacitancia en Farads (100 nanoFarads)
h = 1e-9;      // Paso de tiempo (Dt) en segundos
t_simul = 4e-3; // Tiempo de simulación en segundos (4 ms)

// --- 2. Inicialización de variables ---
t_vector = 0:h:t_simul; // Vector de tiempo. Scilab ya incluye el punto final.
num_points = length(t_vector); // Obtener el número de puntos

// Vectores para almacenar los resultados (inicializados con ceros por defecto)
x1_current = zeros(1, num_points);    // Corriente (i_L)
x2_voltage = zeros(1, num_points);    // Voltaje en capacitor (v_c)
u_input = zeros(1, num_points);       // Voltaje de entrada (V_a)

// Estado inicial [i_L(0); v_c(0)]
state_vector = [0; 0]; // En Scilab, es un vector columna por defecto

// --- 3. Simulación con Euler ---
Va_current = 12; // Voltaje inicial (empezamos con +12V)
for ii = 1:num_points
  current_t = (ii - 1) * h; // Scilab usa índices base 1, como R. (ii-1)*h para obtener el tiempo '0' en el primer paso

  // Cambio de polaridad cada 1 ms
  if (current_t >= 1e-3 && current_t < 2e-3) then
    Va_current = -12;
  elseif (current_t >= 2e-3 && current_t < 3e-3) then
    Va_current = 12;
  elseif (current_t >= 3e-3 && current_t < 4e-3) then
    Va_current = -12;
  end

  // Almacenamiento de resultados en el paso actual antes de la actualización
  t_vector(ii) = current_t;
  u_input(ii) = Va_current;
  x1_current(ii) = state_vector(1);   // Corriente i_L
  x2_voltage(ii) = state_vector(2);   // Voltaje v_c

  // Ecuaciones diferenciales (dx/dt = A*x + B*u)
  di_dt = (-R/L) * state_vector(1) - (1/L) * state_vector(2) + (1/L) * Va_current;
  dvc_dt = (1/C) * state_vector(1);

  // Integración del sistema (Método de Euler)
  state_vector(1) = state_vector(1) + di_dt * h;
  state_vector(2) = state_vector(2) + dvc_dt * h;
end

// --- 4. Gráficos en una sola ventana (TRES Subplots) ---
// Crea o selecciona la figura número 1 para todos los gráficos
hf_main = scf(1); // Selecciona la figura 1
clf();            // Limpia la figura
disp("Creando los tres gráficos en una sola ventana...");
hf_main.background = -2; // Fondo de la figura (ventana) a blanco

// Definir las posiciones de las líneas verticales para los cambios de polaridad
// En Scilab, plot2d permite dibujar líneas directamente en el gráfico
vertical_lines_x = [1e-3, 2e-3, 3e-3];
vertical_lines_y = [-max(abs(x1_current)), max(abs(x1_current))]; // Rango Y para las líneas verticales

// --- Subplot 1: Corriente en el Inductor (i_L(t)) ---
subplot(3, 1, 1); // Activa el panel superior (3 filas, 1 columna, panel 1)
gca().background = -2; // Fondo de los ejes a blanco
plot(t_vector, x1_current * 1000, 'b-', 'LineWidth', 0.8); // Corriente a mA

// Usamos un vector de cadenas para el título multi-línea
title(["Sistema Sobreamortiguado - Entrada escalón variable +-12 V";           "Corriente en el Inductor (mA)"]);
//ylabel('$i_L$ (mA)'); // LaTeX para i_L
xgrid(1); // Activa la cuadrícula
// Añadir líneas verticales para los cambios de polaridad (usando plot2d para dibujo directo)
//for i = 1:length(vertical_lines_x)
  // plot2d(vertical_lines_x(i) * ones(1,2), vertical_lines_y_u, style=2, strf="000", axesflag=0); // Línea discontinua 
//end

// --- Subplot 2: Voltaje en el Capacitor (V_C(t)) ---
subplot(3, 1, 2); // Activa el panel del medio (3 filas, 1 columna, panel 2)
gca().background = -2; // Fondo de los ejes a blanco
plot(t_vector, x2_voltage, 'r-', 'LineWidth', 0.8);
title('Voltaje en el Capacitor (V)'); // LaTeX para v_c
//ylabel('$v_c$ (V)'); // LaTeX para v_c
xgrid(1); // Activa la cuadrícula
// Añadir líneas verticales para los cambios de polaridad
vertical_lines_y_vc = [-max(abs(x2_voltage)), max(abs(x2_voltage))];
//for i = 1:length(vertical_lines_x)
 //   plot2d(vertical_lines_x(i) * ones(1,2), vertical_lines_y_u, style=2, strf="000", axesflag=0); // Línea discontinua
//end

// --- Subplot 3: Entrada Va(t) ---
subplot(3, 1, 3); // Activa el panel inferior (3 filas, 1 columna, panel 3)
gca().background = -2; // Fondo de los ejes a blanco
plot(t_vector, u_input, 'g-', 'LineWidth', 0.8);
title('Voltaje de Entrada (V)'); // LaTeX para V_a
xlabel('tiempo (s)'); // La etiqueta del eje X solo en el gráfico inferior
//ylabel('$V_a$ (V)'); // LaTeX para V_a
xgrid(1); // Activa la cuadrícula
// Añadir líneas verticales para los cambios de polaridad
//vertical_lines_y_u = [-max(abs(u_input)), max(abs(u_input))];
//for i = 1:length(vertical_lines_x)
//  plot2d(vertical_lines_x(i) * ones(1,2), vertical_lines_y_u, style=2, strf="000", axesflag=0); // Línea discontinua sin círculos
//end
disp("Simulación y los tres gráficos en una sola ventana completados.");

// *********** RECORDATORIO PARA GUARDAR EL GRÁFICO MANUALMENTE ***********
// Después de ejecutar el script, la ventana del gráfico permanecerá abierta.
// Para guardar como PNG:
// 1. En la ventana del gráfico, ve a Archivo > Exportar a...
// 2. Navega a la carpeta resultados_graficos/ en tu repositorio local.
// 3. Escribe el nombre del archivo: "RLC_Actividad_3.png" (o el que desees para esta actividad).
// 4. Asegúrate de seleccionar "Portable Network Graphics (*.png)" como tipo de archivo.
// 5. Haz clic en Guardar.
