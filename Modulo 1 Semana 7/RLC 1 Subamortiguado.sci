// Curso HERRAMIENTAS DE CÁLCULO Y SIMULACIÓN DE PROCESOS DINÁMICOS
// Silvia Miró
// Modulo 1 - Semana 7
// Actividad 1: Circuito RLC Subamortiguado

// --- 1. Parámetros del circuito ---
R = 1;      // Resistencia en Ohms
L = 1;      // Inductancia en Henrys
C = 1;      // Capacitancia en Faradays
h = 1e-2;   // Paso de tiempo (Dt)
t_simul = 3; // Tiempo de simulación en segundos.

// --- 2. Matrices del sistema en espacio de estados ---
// Las variables de estado son:
// x(1) = i_L(t)  (Corriente en el inductor)
// x(2) = V_C(t)  (Voltaje en el capacitor)
// La ecuación diferencial es: dx/dt = A*x + B*u
A = [-R/L, -1/L;
     1/C,  0];
B = [1/L;
     0];

// --- 3. Inicialización para la simulación con Euler ---
t = 0:h:t_simul;
x1_euler = zeros(1, length(t));
x2_euler = zeros(1, length(t));
u_entrada = zeros(1, length(t));

// Condiciones iniciales del estado x = [i_L(0); V_C(0)]
x = [0; 0];
Va = 12;

// --- 4. Simulación con Método de Euler ---
for ii = 1:length(t)
    xp = A * x + B * Va;
    x = x + xp * h;
    x1_euler(ii) = x(1);
    x2_euler(ii) = x(2);
    u_entrada(ii) = Va;
end

// --- 5. Solución analítica (para un sistema subamortiguado) ---
omega_0 = 1 / sqrt(L * C);
zeta = R / (2 * sqrt(L / C));
omega_d = omega_0 * sqrt(1 - zeta^2);
i_analitico = Va / (L * omega_d) * exp(-zeta * omega_0 * t) .* sin(omega_d * t);
Vc_analitico = Va * (1 - (omega_0 / omega_d) * exp(-zeta * omega_0 * t) .* sin(omega_d * t + %pi/3));

// --- 6. Gráficos en una sola ventana (Subplots) ---
// Figura número 1 para todos los gráficos del RLC
hf_rlc = scf(1); // Esta línea abre o selecciona la figura
hf_rlc.visible = 'on'; // Aseguramos que la ventana sea visible
clf(); // Limpia la figura
hf_rlc.background = -2; // Fondo de la figura (ventana) a blanco

// --- Subplot 1: Corriente en el Inductor (i_L(t)) ---
subplot(2, 1, 1); // ACTIVA el panel 1 (superior) de una cuadrícula de 2x1
gca().background = -2; // Fondo de los ejes a blanco
plot(t, x1_euler, 'b-', 'LineWidth', 1.5);
plot(t, i_analitico, 'r--', 'LineWidth', 1.5);
// ********************** MODIFICACIÓN AQUÍ **********************
// Usamos un vector de cadenas para el título multi-línea
title(["Sistema Subamortiguado"; "Corriente en el Inductor (A)"]);
// ***************************************************************
//ylabel('Corriente (A)');
legend(['Simulación (Euler)'; 'Solución analítica'], 4);
xgrid(1); // xgrid(1) para activar la cuadrícula

// --- Subplot 2: Voltaje en el Capacitor (V_C(t)) ---
subplot(2, 1, 2); //panel 2 (inferior) de una cuadrícula de 2x1
gca().background = -2; // Fondo de los ejes a blanco
plot(t, x2_euler, 'b-', 'LineWidth', 1.5);
plot(t, Vc_analitico, 'r--', 'LineWidth', 1.5);
title('Voltaje en el Capacitor (V)'); // Título del subplot individual
xlabel('tiempo (s)');
//ylabel('Voltaje (V)');
legend(['Simulación (Euler)'; 'Solución analítica'], 4);
xgrid(1); // xgrid(1) para activar la cuadrícula
disp("Simulación y comparación completadas");

// *********** RECORDATORIO PARA GUARDAR EL GRÁFICO MANUALMENTE ***********
// Después de ejecutar el script, la ventana del gráfico permanecerá abierta.
// Para guardar como PNG:
// 1. En la ventana del gráfico, ve a Archivo > Exportar a...
// 2. Navega a la carpeta resultados_graficos/ en tu repositorio local.
// 3. Escribe el nombre del archivo: "RLC_Actividad_1.png" (o el que desees para esta actividad).
// 4. Asegúrate de seleccionar "Portable Network Graphics (*.png)" como tipo de archivo.
// 5. Haz clic en Guardar.
