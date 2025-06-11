// Curso HERRAMIENTAS DE CÁLCULO Y SIMULACIÓN DE PROCESOS DINÁMICOS
// Silvia Miró
// Modulo 1 - Semana 1
// Actividad 1

// --- 1. Parámetros del circuito ---
R = 1;      // Resistencia en Ohms
L = 1;      // Inductancia en Henrys
C = 1;      // Capacitancia en Farads
h = 1e-2;   // Paso de tiempo (Dt)
t_simul = 3; // Tiempo de simulación en segundos. Ajustado a 8s para mejor visualización.

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

// Crea o selecciona la figura número 1 para todos los gráficos
hf = scf(1);
clf();
disp("Creando gráficos en una sola ventana...");
hf.background = -2; // Fondo de la figura (ventana) a blanco

// --- Subplot 1: Corriente en el Inductor (i_L(t)) ---
subplot(2, 1, 1); // ACTIVA el panel 1 (superior) de una cuadrícula de 2x1
gca().background = -2; // Fondo de los ejes a blanco
plot(t, x1_euler, 'b-', 'LineWidth', 1.5);
plot(t, i_analitico, 'r--', 'LineWidth', 1.5);
title('Corriente en el Inductor (i_L(t))');
ylabel('Corriente (A)');
legend(['Simulación (Euler)'; 'Solución analítica'], 4);
xgrid(1); // <--- CAMBIO AQUÍ: Usamos xgrid(1) para activar la cuadrícula

// --- Subplot 2: Voltaje en el Capacitor (V_C(t)) ---
subplot(2, 1, 2); // ACTIVA el panel 2 (inferior) de una cuadrícula de 2x1
gca().background = -2; // Fondo de los ejes a blanco
plot(t, x2_euler, 'b-', 'LineWidth', 1.5);
plot(t, Vc_analitico, 'r--', 'LineWidth', 1.5);
title('Voltaje en el Capacitor (V_C(t))');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
legend(['Simulación (Euler)'; 'Solución analítica'], 4);
xgrid(1); // <--- CAMBIO AQUÍ: Usamos xgrid(1) para activar la cuadrícula

disp("Simulación y comparación completadas. ¡Deberías ver ambos gráficos en una única ventana!");
