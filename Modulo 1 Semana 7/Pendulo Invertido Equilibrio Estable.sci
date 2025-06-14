// Script en Scilab para la simulación del modelo de Péndulo Invertido sobre Carro.
// Simula el sistema no lineal y su aproximación linealizada.
// ***Equilibrio Estable***
// Modulo 1 Semana 7
// Actividades 1, 2, 3 y 5 para el equilibrio estable

// --- 1. Definir los valores numéricos de los parámetros ---
// m: Masa del péndulo (kg)
// long: Longitud del péndulo (m)

//m = 0.1; long = 0.6;               // Actividad 1 - Eq Estable
//m = 0.2; long = 0.6             // Actividad 2 - Eq Estable
m = 0.01; long = 1.2           // Actividad 3 - Eq Estable

Fricc = 0.1;// Fricc: Fricción del carro (N*s/m) 
g = 9.8;    // Aceleración de la gravedad (m/s^2)
M = 0.5;    // Masa del carro (kg)

mprintf("Parámetros asignados:\n"); 
mprintf("m = %g, Fricc = %g, long = %g, g = %g, M = %g\n", m, Fricc, long, g, M);

// --- 2. Parámetros de la simulación ---
h = 0.0001; // Paso de tiempo (s)
tiempo_total = 10; // Tiempo total de simulación en segundos
num_pasos = int(tiempo_total / h); // Número total de pasos
t = 0:h:(num_pasos*h); // Vector de tiempo
n_steps = length(t);

mprintf("Paso de tiempo (h): %g s, Tiempo total: %g s, Número de pasos: %d\n", h, tiempo_total, n_steps);

// --- 3. Inicialización de arrays para la simulación NO lineal ---
omega = zeros(1, n_steps);
alfa = zeros(1, n_steps);
p = zeros(1, n_steps);
p_p = zeros(1, n_steps);
u = zeros(1, n_steps);

// --- 4. Condición inicial para la simulación NO lineal ---
alfa(1) = 3.2; // Ángulo inicial 
omega(1) = 0;
p(1) = 0;
p_p(1) = 0;

mprintf("\nCondiciones iniciales NO lineales: alfa = %g, omega = %g, p = %g, p_p = %g\n", alfa(1), omega(1), p(1), p_p(1));

// --- 5. Matrices del sistema LINEALIZADO (Sontag, equilibrio estable) ---
Mat_A = [
    0, 1, 0, 0;
    0, -Fricc/M, -m*g/M, 0;
    0, 0, 0, 1;
    0, -Fricc/(long*M), -g*(m+M)/(long*M), 0
];

Mat_B = [0; 1/M; 0; 1/(long*M)];

mprintf("\nMatriz A (linealizada):\n"); // MODIFICADO
disp(Mat_A); // disp para matrices es el más adecuado
mprintf("Matriz B (linealizada):\n"); // MODIFICADO
disp(Mat_B); // disp para matrices es el más adecuado

// --- 6. Punto de equilibrio para la linealización y condición inicial para sim. lineal ---
X0_equilibrium = [0; 0; %pi; 0];

// Condición inicial para el estado lineal (x_lineal).
x_lineal = [0; 0; alfa(1); 0];

// --- 7. Inicialización de arrays para la simulación LINEALIZADA ---
pl = zeros(1, n_steps);
p_pl = zeros(1, n_steps);
alfal = zeros(1, n_steps);
omegal = zeros(1, n_steps);

// --- 8. Simulación (bucle principal) ---
mprintf("\nIniciando simulación...\n"); // MODIFICADO

for i = 1:(n_steps - 1)
    // --- Simulación NO LINEAL ---
    sin_alfa = sin(alfa(i));
    cos_alfa = cos(alfa(i));

    // Matriz de coeficientes para el sistema acoplado [p_pp, alfa_pp]
    A_no_lin = [
        1, m*long*cos_alfa/(M+m);
        cos_alfa/long, 1
    ];

    // Vector de términos independientes
    b_no_lin = [
        (u(i) + m*long*omega(i)^2*sin_alfa - Fricc*p_p(i))/(M+m);
        g*sin_alfa/long
    ];

    // Resolver el sistema lineal para las aceleraciones
    sol_no_lin = A_no_lin \ b_no_lin;
    p_pp_val = sol_no_lin(1);
    alfa_pp_val = sol_no_lin(2);

    // Euler NO lineal: Actualizar estados
    p_p(i+1) = p_p(i) + h * p_pp_val;
    p(i+1) = p(i) + h * p_p(i);
    omega(i+1) = omega(i) + h * alfa_pp_val;
    alfa(i+1) = alfa(i) + h * omega(i);

    // --- Simulación LINEALIZADA ---
    xp_lineal = Mat_A * (x_lineal - X0_equilibrium) + Mat_B * u(i);

    // Euler lineal: Actualizar estado lineal
    x_lineal = x_lineal + h * xp_lineal;

    // Guardar estados LINEALIZADOS para graficar
    pl(i) = x_lineal(1);
    p_pl(i) = x_lineal(2);
    alfal(i) = x_lineal(3);
    omegal(i) = x_lineal(4);
end

// Guardar el último paso para los arrays de resultados linealizados
pl(n_steps) = x_lineal(1);
p_pl(n_steps) = x_lineal(2);
alfal(n_steps) = x_lineal(3);
omegal(n_steps) = x_lineal(4);

// --- 9. Gráficos ---
mprintf("\nGenerando gráficos...\n"); 

// Figura 1: Comportamiento de los estados y acción de control
hf_main1 = scf(1);
clf();
hf_main1.background = -2;

subplot(3, 2, 1);
gca().background = -2;
plot(t, omega, 'b-', 'LineWidth', 1);
plot(t, omegal, 'k--', 'LineWidth', 1);
title("Velocidad Angular (rad/s)");
xlabel("tiempo (s)")
//ylabel("Vel. Angular (rad/s)");
xgrid(1); legend(["No lineal", "Lineal"]);

subplot(3, 2, 2);
gca().background = -2;
plot(t, alfa, 'b-', 'LineWidth', 1);
plot(t, %pi * ones(1, n_steps), 'k:', 'LineWidth', 0.8);
plot(t, alfal, 'k--', 'LineWidth', 1);
title("Ángulo (rad)");
xlabel("tiempo (s)")
//ylabel("Ángulo (rad)");
xgrid(1); legend(["No lineal", "Lineal"]);

subplot(3, 2, 3);
gca().background = -2;
plot(t, p, 'b-', 'LineWidth', 1);
plot(t, pl, 'k--', 'LineWidth', 1);
title("Posición Carro (m)");
xlabel("tiempo (s)")
//ylabel("Posición (m)");
xgrid(1); legend(["No lineal", "Lineal"]);

subplot(3, 2, 4);
gca().background = -2;
plot(t, p_p, 'b-', 'LineWidth', 1);
plot(t, p_pl, 'k--', 'LineWidth', 1);
title("Velocidad Carro (m/s)");
xlabel("tiempo (s)")
//ylabel("Velocidad (m/s)");
xgrid(1); legend(["No lineal", "Lineal"]);

subplot(3, 1, 3);
gca().background = -2;
plot(t, u, 'b-', 'LineWidth', 1);
title("Acción de Control u (N)");
xlabel("tiempo (s)");
//ylabel("Fuerza (N)");
xgrid(1);


// Figura 2: Diagramas de fase
hf_main2 = scf(2);
clf();
hf_main2.background = -2;

subplot(2, 2, 1);
gca().background = -2;
plot(alfa, omega, 'b-', 'LineWidth', 1);
plot(alfal, omegal, 'k--', 'LineWidth', 1);
title("Espacio de Fase: Ángulo vs. Vel. Angular");
xlabel("Ángulo (rad)");
ylabel("Vel. Angular (rad/s)");
xgrid(1); legend(["No lineal", "Lineal"]);

subplot(2, 2, 2);
gca().background = -2;
plot(p, p_p, 'b-', 'LineWidth', 1);
plot(pl, p_pl, 'k--', 'LineWidth', 1);
title("Espacio de Fase: Posición vs. Vel. Carro");
xlabel("Posición (m)");
ylabel("Velocidad (m/s)");
xgrid(1); legend(["No lineal", "Lineal"]);

// Comentado para no guardar los datos en un archivo .mat
// save('Datos_Controlador.mat', 'alfa', 'omega', 'p', 'p_p', 'u', 'pl', 'p_pl', 'alfal', 'omegal', 't');

// Línea final de disp 
mprintf('Simulación completada y datos listos para visualizar.\n'); 
