// Curso HERRAMIENTAS DE CÁLCULO Y SIMULACIÓN DE PROCESOS DINÁMICOS
// Silvia Miró
// Modulo 1 - Semana 7
// Actividad 2: Sitema RLC Sobreamortiguado

// --- 1. Parámetros del circuito ---
R = 2.2e3;     // Resistencia en Ohms (2.2 kOhms)
L = 10e-6;     // Inductancia en Henrys (10 microHenrys)
C = 100e-9;    // Capacitancia en Farads (100 nanoFarads)
h = 1e-9;      // Paso de tiempo (Dt) en segundos
t_simul = 1e-3; // Tiempo de simulación en segundos (1 milisegundo)

// --- 2. Inicialización ---
t = 0:h:t_simul; // Vector de tiempo
// Inicialización de vectores para almacenar resultados
x1_euler = zeros(1, length(t)); // Corriente i_L(t)
x2_euler = zeros(1, length(t)); // Voltaje V_C(t)
u_entrada = zeros(1, length(t)); // Entrada Va(t)

// Condiciones iniciales
x = [0; 0]; // Vector de estado [i_L(0); V_C(0)]
Va = 12;    // Valor del escalón de tensión de 12V

// --- 3. Simulación con Método de Euler ---
for ii = 1:length(t)
    // Calcula las derivadas en el punto actual (di_dt, dvc_dt)
    di_dt = (-R/L) * x(1) - (1/L) * x(2) + (1/L) * Va; // Ecuación para di_L/dt
    dvc_dt = (1/C) * x(1); // Ecuación para dV_C/dt

    // Método de Euler: x(k+1) = x(k) + dx/dt(k) * h
    x(1) = x(1) + di_dt * h; // Actualiza la corriente del inductor
    x(2) = x(2) + dvc_dt * h; // Actualiza el voltaje en el capacitor

    // Almacena los resultados para graficar
    x1_euler(ii) = x(1); // Almacena corriente
    x2_euler(ii) = x(2); // Almacena voltaje en C
    u_entrada(ii) = Va;  // Almacena entrada
end

// --- 4. Solución analítica (PARA UN SISTEMA SOBREAMORTIGUADO) ---
// NOTA IMPORTANTE: Con estos valores R, L, C, el circuito es **sobreamortiguado**.
// Hemos corregido las fórmulas analíticas para que coincidan con este caso.
alpha = R / (2 * L);             // Coeficiente de amortiguamiento
omega_0 = 1 / sqrt(L * C);       // Frecuencia natural no amortiguada

// Cálculo de las raíces para el caso sobreamortiguado
// Usamos sqrt(alpha^2 - omega_0^2) ya que alpha > omega_0
s1_val = -alpha + sqrt(alpha^2 - omega_0^2);
s2_val = -alpha - sqrt(alpha^2 - omega_0^2);

// Solución analítica para el voltaje en el capacitor V_C(t)
Vc_analitico = Va * (1 - (s1_val / (s1_val - s2_val)) * exp(s2_val * t) + (s2_val / (s1_val - s2_val)) * exp(s1_val * t));

// Solución analítica para la corriente en el inductor i_L(t)
i_analitico = (Va * C * s1_val * s2_val / (s1_val - s2_val)) * (exp(s1_val * t) - exp(s2_val * t));

// --- 5. Gráficos en una sola ventana (TRES Subplots) ---
// Crea o selecciona la figura número 1 para todos los gráficos
hf_main = scf(1); // Selecciona la figura 1
clf();            // Limpia la figura
disp("Creando los tres gráficos en una sola ventana...");
hf_main.background = -2; // Fondo de la figura (ventana) a blanco

// --- Subplot 1: Corriente en el Inductor (i_L(t)) ---
subplot(3, 1, 1); // **ACTIVA el panel superior (3 filas, 1 columna, panel 1)**
gca().background = -2; // Fondo de los ejes a blanco
// No necesitamos real() porque los resultados analíticos ahora deberían ser reales
plot(t, x1_euler, 'b-', 'LineWidth', 1.5);
plot(t, i_analitico, 'r--', 'LineWidth', 1.5);
// Usamos un vector de cadenas para el título multi-línea
title(["Sistema Sobreamortiguado";       "Corriente en el Inductor (A)"]);
//ylabel('Corriente (A)');
legend(['Simulación (Euler)'; 'Solución analítica'], 4);
xgrid(1); // Activa la cuadrícula

// --- Subplot 2: Voltaje en el Capacitor (V_C(t)) ---
subplot(3, 1, 2); // **ACTIVA el panel del medio (3 filas, 1 columna, panel 2)**
gca().background = -2; // Fondo de los ejes a blanco
// No necesitamos real() porque los resultados analíticos ahora deberían ser reales
plot(t, x2_euler, 'b-', 'LineWidth', 1.5);
plot(t, Vc_analitico, 'r--', 'LineWidth', 1.5);
title('Voltaje en el Capacitor (V)'); // Título del subplot individual
//ylabel('Voltaje (V)');
legend(['Simulación (Euler)'; 'Solución analítica'], 4);
xgrid(1); // Activa la cuadrícula

// --- Subplot 3: Entrada Va(t) ---
subplot(3, 1, 3); // **ACTIVA el panel inferior (3 filas, 1 columna, panel 3)**
gca().background = -2; // Fondo de los ejes a blanco
plot(t, u_entrada, 'g-', 'LineWidth', 1.5);
title('Entrada de Tensión (V)'); // Título del subplot individual
xlabel('tiempo (s)'); // La etiqueta del eje X solo en el gráfico inferior
//ylabel('Voltaje (V)');
xgrid(1); // Activa la cuadrícula
legend('Entrada Va(V)', 4);
disp("Simulación y los tres gráficos en una sola ventana completados.");

// *********** RECORDATORIO PARA GUARDAR EL GRÁFICO MANUALMENTE ***********
// Después de ejecutar el script, la ventana del gráfico permanecerá abierta.
// Para guardar como PNG:
// 1. En la ventana del gráfico, ve a Archivo > Exportar a...
// 2. Navega a la carpeta resultados_graficos/ en tu repositorio local.
// 3. Escribe el nombre del archivo: "RLC_Actividad_2.png" (o el que desees para esta actividad).
// 4. Asegúrate de seleccionar "Portable Network Graphics (*.png)" como tipo de archivo.
// 5. Haz clic en Guardar.
