// Script en Scilab para la simulación del modelo de Péndulo Invertido sobre Carro.
// Silvia Miro
// Modulo 1 - Semana 7
// Pendulo Invertido Equilibrio Inestable
// Actividades 1, 2, 3 y 5

// --- 1. Definir los valores numéricos de los parámetros ---
 // m_val: Masa del péndulo (kg), long_val: Longitud del péndulo (m)

//m_val = 0.1; long_val = 0.6;     //  Actividad 1
//m_val = 0.2; long_val = 0.6;     //  Actividad 2
m_val = 0.01; long_val = 1.2;     //  Actividad 3

Fricc_val = 0.1;   // Coeficiente de fricción
g_val = 9.8;     // Aceleración de la gravedad (m/s^2)
M_val = 0.5;     // Masa del carro (kg)
h_val = 0.0001;  // Paso de tiempo (s)
tiempo_final = 10; // Tiempo total de simulación (s)
n_steps = int(tiempo_final / h_val); // Número de pasos

mprintf("--- Parámetros del Sistema ---\n");
mprintf("Masa Péndulo (m): %.2f kg\n", m_val);
mprintf("Fricción Carro (Fricc): %.2f N*s/m\n", Fricc_val);
mprintf("Longitud Péndulo (l): %.2f m\n", long_val);
mprintf("Gravedad (g): %.1f m/s^2\n", g_val);
mprintf("Masa Carro (M): %.2f kg\n", M_val);
mprintf("Paso de tiempo (h): %g s\n", h_val);
mprintf("Tiempo total simulación: %d s\n", tiempo_final);
mprintf("Número de pasos: %d\n", n_steps);
mprintf("----------------------------\n\n");

// --- 2. Inicialización de arrays ---
omega = zeros(1, n_steps + 1);     // Velocidad angular del péndulo
alfa = zeros(1, n_steps + 1);      // Ángulo del péndulo
p = zeros(1, n_steps + 1);         // Posición del carro
p_p = zeros(1, n_steps + 1);       // Velocidad del carro
u = zeros(1, n_steps + 1);         // Acción de control (entrada), inicializada a ceros por defecto

// --- 3. Condiciones iniciales (equilibrio inestable alfa ≈ 0) ---
alfa(1) = -0.01; // Pequeño desvío desde alfa = 0 (vertical hacia arriba)
omega(1) = 0;
p(1) = 0;
p_p(1) = 0;

mprintf("--- Condiciones Iniciales ---\n");
mprintf("Ángulo inicial (alfa(1)): %.4f rad\n", alfa(1));
mprintf("Velocidad angular inicial (omega(1)): %.4f rad/s\n", omega(1));
mprintf("Posición inicial carro (p(1)): %.4f m\n", p(1));
mprintf("Velocidad inicial carro (p_p(1)): %.4f m/s\n", p_p(1));
mprintf("----------------------------\n\n");

// --- 4. Matrices del sistema linealizado (equilibrio inestable alfa ≈ 0) ---
// La matriz A está linealizada alrededor de alfa = 0 (péndulo vertical hacia arriba).
// El estado es [p, p_dot, alfa, alfa_dot]
Mat_A = [
    0, 1, 0, 0;
    0, -Fricc_val / M_val, -m_val * g_val / M_val, 0;
    0, 0, 0, 1;
    0, Fricc_val / (long_val * M_val), (M_val + m_val) * g_val / (long_val * M_val), 0
];

Mat_B = [0; 1 / M_val; 0; 1 / (long_val * M_val)];

mprintf("--- Matrices del Sistema Linealizado (alfa ≈ 0) ---\n");
disp(Mat_A);
mprintf("\n");
disp(Mat_B);
mprintf("----------------------------------------------\n\n");

// --- 5. Punto de equilibrio para la linealización ---
X0_equilibrium = [0; 0; 0; 0]; // [posición_carro, velocidad_carro, ángulo_pendulo, vel_angular_pendulo]
                               // Equilibrio: carro en 0, quieto; péndulo en 0 (vertical hacia arriba), quieto.

// --- 6. Inicialización resultados del sistema linealizado ---
x_linear = [p(1); p_p(1); alfa(1); omega(1)]; // Estado inicial del sistema linealizado (columna vector)
pl = zeros(1, n_steps + 1);    // Posición carro linealizada
p_pl = zeros(1, n_steps + 1);  // Velocidad carro linealizada
alfal = zeros(1, n_steps + 1); // Ángulo linealizado
omegal = zeros(1, n_steps + 1); // Velocidad angular linealizada

// --- 7. Simulación ---
mprintf("Iniciando simulación...\n");
tic(); // Iniciar el temporizador

for i = 1:n_steps
    // --- Simulación No Lineal ---
    sin_phi = sin(alfa(i));
    cos_phi = cos(alfa(i));

    A_nl = [
        1, m_val * long_val * cos_phi / (M_val + m_val);
        cos_phi / long_val, 1
    ];

    b_nl = [
        (u(i) + m_val * long_val * omega(i)^2 * sin_phi - Fricc_val * p_p(i)) / (M_val + m_val);
        g_val * sin_phi / long_val
    ];

    // Resolver el sistema lineal para las aceleraciones
    sol_nl = A_nl \ b_nl;

    p_pp_val = sol_nl(1);
    tita_pp_val = sol_nl(2); // Aceleración angular

    // Integración de Euler para el sistema No Lineal
    p_p(i + 1) = p_p(i) + h_val * p_pp_val;
    p(i + 1) = p(i) + h_val * p_p(i);
    omega(i + 1) = omega(i) + h_val * tita_pp_val;
    alfa(i + 1) = alfa(i) + h_val * omega(i);

    // --- Simulación Lineal ---
    // Guardar estados del paso actual antes de la actualización
    pl(i) = x_linear(1);
    p_pl(i) = x_linear(2);
    alfal(i) = x_linear(3);
    omegal(i) = x_linear(4);

    xp_linear = Mat_A * (x_linear - X0_equilibrium) + Mat_B * u(i);
    x_linear = x_linear + h_val * xp_linear;
end

// Asignar los últimos valores calculados del sistema linealizado
pl(n_steps + 1) = x_linear(1);
p_pl(n_steps + 1) = x_linear(2);
alfal(n_steps + 1) = x_linear(3);
omegal(n_steps + 1) = x_linear(4);

// --- 8. Vector de tiempo para los gráficos ---
t_vector = 0:h_val:(n_steps*h_val);

// --- 9. Gráficos ---
mprintf("Generando gráficos...\n");

// --- Submuestreo para Gráficos ---
max_points_for_plotting = 1000; // Máximo de puntos a mostrar en los gráficos por línea
plot_indices = 1:(n_steps + 1); // Por defecto, todos los puntos

if (n_steps + 1) > max_points_for_plotting then
    step_for_plotting = ceil((n_steps + 1) / max_points_for_plotting);
    plot_indices = 1:step_for_plotting:(n_steps + 1);
    mprintf("Submuestreando datos para gráficos. Se mostrarán aproximadamente %d puntos por línea.\n", length(plot_indices));
end

// Variables subsampleadas para graficar
t_plot = t_vector(plot_indices);
omega_plot = omega(plot_indices);
alfa_plot = alfa(plot_indices);
p_plot = p(plot_indices);
p_p_plot = p_p(plot_indices);
u_plot = u(plot_indices);

omegal_plot = omegal(plot_indices);
alfal_plot = alfal(plot_indices);
pl_plot = pl(plot_indices);
p_pl_plot = p_pl(plot_indices);


// Definir límites de Y específicos
y_lim_alfa_min = -6.0;
y_lim_alfa_max = 6.0;
y_lim_omega_min = -10.0;
y_lim_omega_max = 10.0;
y_lim_p_min = -0.2;
y_lim_p_max = 0.2;
y_lim_pp_min = -1.5;
y_lim_pp_max = 1.5;

// --- Figura 1: Comportamiento de los estados y acción de control (VISIBLE) ---
hf_main1 = scf(1); // Crea/selecciona la figura
hf_main1.visible = 'on'; 
clf();
hf_main1.background = -2; // Fondo blanco

subplot(3, 2, 1);
gca().background = -2;
plot(t_plot, omega_plot, 'b-', 'LineWidth', 1);
plot(t_plot, omegal_plot, 'k--', 'LineWidth', 1);
title("Velocidad Angular (rad/s)");
xlabel("tiempo (s)")
//ylabel("Velocidad Angular (rad/s)");
xgrid(1); legend(["No lineal", "Lineal"]);
gca().data_bounds = [min(t_plot), y_lim_omega_min; max(t_plot), y_lim_omega_max]; // Set Y limits

subplot(3, 2, 2);
gca().background = -2;
plot(t_plot, alfa_plot, 'b-', 'LineWidth', 1);
plot(t_plot, zeros(1, length(t_plot)), 'k:', 'LineWidth', 0.8); // Línea de equilibrio en 0 (longitud ajustada)
plot(t_plot, %pi * ones(1, length(t_plot)), 'g:', 'LineWidth', 0.8); // Línea de equilibrio en pi (longitud ajustada)
plot(t_plot, -%pi * ones(1, length(t_plot)), 'g:', 'LineWidth', 0.8); // Línea de equilibrio en -pi (longitud ajustada)
plot(t_plot, alfal_plot, 'k--', 'LineWidth', 1);
title("Ángulo (rad)");
xlabel("tiempo (s)")
//ylabel("Ángulo (rad)");
xgrid(1); legend(["No lineal", "Lineal"]);
gca().data_bounds = [min(t_plot), y_lim_alfa_min; max(t_plot), y_lim_alfa_max]; // Set Y limits

subplot(3, 2, 3);
gca().background = -2;
plot(t_plot, p_plot, 'b-', 'LineWidth', 1);
plot(t_plot, pl_plot, 'k--', 'LineWidth', 1);
title("Posición carro (m)");
xlabel("tiempo (s)")
//ylabel("Posición (m)");
xgrid(1); legend(["No lineal", "Lineal"]);
gca().data_bounds = [min(t_plot), y_lim_p_min; max(t_plot), y_lim_p_max]; // Set Y limits

subplot(3, 2, 4);
gca().background = -2;
plot(t_plot, p_p_plot, 'b-', 'LineWidth', 1);
plot(t_plot, p_pl_plot, 'k--', 'LineWidth', 1);
title("Velocidad carro (m/s)");
xlabel("tiempo (s)")
//ylabel("Velocidad (m/s)");
xgrid(1); legend(["No lineal", "Lineal"]);
gca().data_bounds = [min(t_plot), y_lim_pp_min; max(t_plot), y_lim_pp_max]; // Set Y limits

subplot(3, 1, 3);
gca().background = -2;
plot(t_plot, u_plot, 'b-', 'LineWidth', 1);
title("Acción de control (N)");
xlabel("tiempo (s)");
//ylabel("Fuerza (N)");
xgrid(1);

// xsave("estados_simulacion.png", hf_main1); // Comentado para que no se guarde automáticamente (si lo quieres, descomenta)
// close(hf_main1); // ¡ELIMINADO para que la figura 1 permanezca abierta!

// --- Figura 2: Gráfico de Fase: Ángulo vs. Vel. Angular (VISIBLE) ---
hf_main2 = scf(2); // Crea/selecciona la figura
hf_main2.visible = 'on'; // 
clf();
hf_main2.background = -2;

subplot(2, 1, 1);
gca().background = -2;
plot(alfa_plot, omega_plot, 'b-', 'LineWidth', 0.8);
plot(alfal_plot, omegal_plot, 'k--', 'LineWidth', 0.8);
title("Espacio de Fases: Ángulo vs. Vel. Angular");
xlabel("Ángulo (rad)");
ylabel("Velocidad angular (rad/s)");
xgrid(1); legend(["No lineal", "Lineal"]);
gca().data_bounds = [y_lim_alfa_min, y_lim_omega_min; y_lim_alfa_max, y_lim_omega_max];

// xsave("pendulo_fases_angulo.png", hf_main2); // Comentado para que no se guarde automáticamente
// close(hf_main2); // 

// --- Figura 3: Gráfico de Fase: Posición vs. Vel. Carro (COMENTADA) ---
// La siguiente sección ha sido comentada para excluir la Figura 3 del proceso.

//hf_main3 = scf(3); // Crea/selecciona la figura
//hf_main3.visible = 'on'; // Hace la figura invisible
//clf();
//hf_main3.background = -2; // Fondo para la imagen guardada

subplot(2, 1, 2); 
gca().background = -2;
plot(p_plot, p_p_plot, 'b.', 'MarkerSize', 2);
plot(pl_plot, p_pl_plot, 'k.', 'MarkerSize', 2);
title("Espacio de Fases: Posición vs. Vel. Carro");
xlabel("Posición carro (m)");
ylabel("Velocidad carro (m/s)");
xgrid(1); legend(["No lineal", "Lineal"], 5);
gca().data_bounds = [y_lim_p_min, y_lim_pp_min; y_lim_p_max, y_lim_pp_max];

//xsave("pendulo_fases_carro.png", hf_main3); // Guarda como PNG
//close(hf_main3); // Cierra la figura para liberar recursos


// --- Tiempo de cómputo ---
end_time_r = toc(); // Obtener el tiempo transcurrido
mprintf("\nTiempo de cálculo = %g segundos\n", end_time_r);

// --- Análisis de Autovalores ---
mprintf("\n--- Análisis de Autovalores para el Equilibrio Inestable ---\n");
// Calcular los autovalores de la matriz A del sistema linealizado
eigenvalues_inestable = spec(Mat_A); // OBTENER SOLO LA PRIMERA SALIDA

// Función para interpretar los autovalores en Scilab
function interpret_eigenvalues_scilab(eigenvalues)
    for i = 1:length(eigenvalues)
        lambda = eigenvalues(i);
        real_part = real(lambda);
        imag_part = imag(lambda);
        mprintf("Autovalor %d: %.4f %+.4fi\n", i, real_part, imag_part);
        if abs(real_part) < 1e-5 & abs(imag_part) < 1e-5 then
            mprintf("  - Comportamiento: Constante (punto de equilibrio o linea)\n");
        elseif abs(real_part) < 1e-5 then
            mprintf("  - Comportamiento: Oscilación Pura (Centro)\n");
            mprintf("  - Frecuencia de Oscilación: %.2f rad/s\n", abs(imag_part));
        elseif real_part < 0 then
            if abs(imag_part) < 1e-5 then
                mprintf("  - Comportamiento: Decaimiento Exponencial (Convergencia Directa)\n");
                mprintf("  - Tasa de Decaimiento: %.2f\n", abs(real_part));
            else
                mprintf("  - Comportamiento: Oscilación Amortiguada (Espiral Convergente)\n");
                mprintf("  - Frecuencia natural de Oscilación (aproximada): %.2f rad/s\n", abs(imag_part));
                mprintf("  - Tasa de Amortiguación: %.2f\n", abs(real_part));
            end
        else // real_part > 0
            if abs(imag_part) < 1e-5 then
                mprintf("  - Comportamiento: Crecimiento Exponencial (Divergencia Directa) - ¡INESTABLE!\n");
                mprintf("  - Tasa de Crecimiento: %.2f\n", real_part);
            else
                mprintf("  - Comportamiento: Oscilación Creciente (Espiral Divergente) - ¡INESTABLE!\n");
                mprintf("  - Frecuencia natural de Oscilación (aproximada): %.2f rad/s\n", abs(imag_part));
                mprintf("  - Tasa de Crecimiento: %.2f\n", real_part);
            end
        end
    end
endfunction

mprintf("\n");
interpret_eigenvalues_scilab(eigenvalues_inestable);

// --- Verificación de Valores ADICIONALES (para ver la divergencia inicial) ---
mprintf("\n--- Verificación de Valores Detallada (Alfa) ---\n");
// Valores a t = 0.00s
mprintf("Valores a t = %.4f s (índice %d):\n", t_vector(1), 1);
mprintf("  Alfa_L (Lineal): %.6f rad\n", alfal(1));
mprintf("  Alfa_NL (No Lineal): %.6f rad\n", alfa(1));

// Valores a t = 0.05s
idx_0_05s = find(abs(t_vector - 0.05) == min(abs(t_vector - 0.05)), 1);
mprintf("Valores a t = %.4f s (índice %d):\n", t_vector(idx_0_05s), idx_0_05s);
mprintf("  Alfa_L (Lineal): %.6f rad\n", alfal(idx_0_05s));
mprintf("  Alfa_NL (No Lineal): %.6f rad\n", alfa(idx_0_05s));

// Valores a t = 0.10s
idx_0_10s = find(abs(t_vector - 0.10) == min(abs(t_vector - 0.10)), 1);
mprintf("Valores a t = %.4f s (índice %d):\n", t_vector(idx_0_10s), idx_0_10s);
mprintf("  Alfa_L (Lineal): %.6f rad\n", alfal(idx_0_10s));
mprintf("  Alfa_NL (No Lineal): %.6f rad\n", alfa(idx_0_10s));

// Valores a t = 0.20s
idx_0_20s = find(abs(t_vector - 0.20) == min(abs(t_vector - 0.20)), 1);
mprintf("Valores a t = %.4f s (índice %d):\n", t_vector(idx_0_20s), idx_0_20s);
mprintf("  Alfa_L (Lineal): %.6f rad\n", alfal(idx_0_20s));
mprintf("  Alfa_NL (No Lineal): %.6f rad\n", alfa(idx_0_20s));
mprintf("----------------------------------------\n");

mprintf("Simulación completada y datos listos para visualizar.\n");
