// Script en Scilab para la simulación del modelo de vuelo longitudinal linealizado
// Silvia Miró
// Modulo 1 - Semana 7
// Actividad 2.2 y 2.3

// --- 1. Parámetros del modelo ---
omega = 0.2;    // Frecuencia natural (rad/s)
a = 0.01;       // Coeficiente de amortiguamiento
b = 2;          // Constante de control
//c = 100; t_final = 5;  // Velocidad de vuelo (m/s), tiempo de simulación (s)- Actividad 2.2
c = 50; t_final = 20;          // Velocidad de vuelo (m/s) - Actividad 2.3

// --- 2. Parámetros de la simulación ---
h_step = 1e-3;              // Paso de tiempo (s)
              
t = 0:h_step:t_final;       // Vector de tiempo
n_steps = length(t);       // Número de pasos

// --- 3. Matrices del sistema lineal ---
// dx/dt = A*x + B*u
// x = [alpha; phi; phi_dot; h]

A = [-a,         a,          0,              0;
     0,          0,          1,              0;
     omega^2, -omega^2,      0,              0;
     c,          0,          0,              0]; 

B = [0;
     0;
     omega^2 * b;
     0];

// --- 4. Autovalores de la matriz A ---
// En Scilab, la función para autovalores es spec()
eigenvalues_A = spec(A);
disp('Autovalores de la matriz A:');
disp(eigenvalues_A);

// --- 5. Condición inicial ---
x0 = [0.1;    // alpha(0) - Angulo de ataque inicial
      0;      // phi(0)   - Angulo de cabeceo inicial
      0;      // phi_dot(0) - Velocidad angular de cabeceo inicial
      0];      // h(0)     - Altitud inicial

// --- 6. Inicialización de la matriz de estados ---
x_results = zeros(4, n_steps);
x_results(:, 1) = x0;

// --- 7. Entrada de control (nula para evaluar el comportamiento inherente) ---
u = zeros(1, n_steps); // Entrada nula como en el script de Octave

// --- 8. Simulación con el método de Euler hacia adelante ---
state_vector_current = x0;

for n = 1:(n_steps - 1)
    // x(n+1) = x(n) + h * (A * x(n) + B * u(n))
    dx_dt = A * state_vector_current + B * u(n);
    state_vector_current = state_vector_current + h_step * dx_dt;
    x_results(:, n + 1) = state_vector_current;
end

// --- 9. Gráficos en una sola ventana (4 Subplots) ---

hf_main = scf(1); // Selecciona la figura 1
clf();            // Limpia la figura
disp("Generando gráficos de la simulación del modelo Octave en Scilab...");
hf_main.background = -2; // Fondo de la figura a blanco

// --- Subplot 1: Ángulo de Ataque (rad) ---
subplot(4, 1, 1);
gca().background = -2;
plot(t, x_results(1,:), 'b-', 'LineWidth', 0.8);
title(msprintf('c=%g m/s: Ángulo de Ataque (rad)', c));
//ylabel('Ángulo (rad)');
xgrid(1);

// --- Subplot 2: Ángulo de Cabeceo (rad) ---
subplot(4, 1, 2);
gca().background = -2;
plot(t, x_results(2,:), 'r-', 'LineWidth', 0.8);
title(msprintf('c=%g m/s: Ángulo de Cabeceo (rad)', c));
//ylabel('Ángulo (rad)');
xgrid(1);

// --- Subplot 3: Velocidad angular de cabeceo (rad/s) ---
subplot(4, 1, 3);
gca().background = -2;
plot(t, x_results(3,:), 'g-', 'LineWidth', 0.8);
title(msprintf('c=%g m/s: Velocidad Angular (rad/s)', c));
//ylabel('Velocidad (rad/s)');
xgrid(1);

// --- Subplot 4: Altitud (m) ---
subplot(4, 1, 4);
gca().background = -2;
plot(t, x_results(4,:), 'm-', 'LineWidth', 0.8);
title(msprintf('c=%g m/s: Altitud (m)', c));
xlabel('tiempo (s)');
//ylabel('Altitud (m)');
xgrid(1);

disp("Simulación del modelo Octave en Scilab completada.");
