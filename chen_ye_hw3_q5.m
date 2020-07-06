% Jiahao (Derek) Ye
% Jiyao (Lucas) Chen
%
% chen_ye_hw3_q4.m
% CS346 -- Computational Modeling & Simulation 
% Spring, 2020


% initialize constants for the simulation
dt = 0.001;                 % ms
simulation_time = 3;        % ms
N = simulation_time/dt;     % number of iterations

V = zeros(1,N);             % potential in mV
t = zeros(1,N);             % time in ms
n = zeros(1,N);             % K activation gating variable
m = zeros(1,N);             % Na activation gating variable
h = zeros(1,N);             % Na inactivation gating variable
K_i = zeros(1,N);             % K+ concentration inside
K_o = zeros(1,N);             % K+ concentration outside
Na_i = zeros(1,N);            % Na+ concentration inside
Na_o = zeros(1,N);            % Na+ concentration outside

% initialize constants for the model
V(1) = -65;          % initial action potential in mV
C = 0.1;             % capacitance in microF/cm3

n(1) = 0.317;        % initial K activation gating variable
m(1) = 0.05;         % initial Na activation gating variable
h(1) = 0.6;          % initial Na inactivation gating
t(1) = 0;            % initial time
K_i(1) = 150;        % K ion initial concentration inside
K_o(1) = 5.5;        % K ion initial concentration outside
Na_i(1) = 15;        % Na ion initial concentration inside
Na_o(1) = 150;       % Na ion initial concentration outside

V_Na_thres = 50;   % threshold for Na_close | K_open in mV
V_Na_start = -55;    % initial condition for Na to open
V_K_disp = -77;      % displacement from the equilibrium for K
V_Na_disp = 50;      % displacement from the equilibrium for Na
V_leak = -54.4;      % displacement for lekage in mV
gK = 36;             % max K conductance in mS/cm^2
gNa = 120;           % max Na conductance in mS/cm^2
gL = 0.3;            % max leakage conductance in mS/cm^2

Iext = 15;          % applied current in nA
stimulus_start = 0.5;     % stimulus start time in ms
stimulus_duration = 0.5;  % stimulus last for this amount of ms

dKdt = -2;          % rate of change of K outside
dNadt = -3;         % rate of change of Na inside

% initialize anonymous functions
% anonymous functions for opening/closing rate constants
an = @(V) 0.01 * (V+55)/(1-exp(-(V+55)/10));
am = @(V) 0.1 * (V+40)/(1-exp(-(V+40)/10));
ah = @(V) 0.07 * exp(-(V+65)/20);
bn = @(V) 0.125 * exp(-(V+65)/80);
bm = @(V) 4 * exp(-(V+65)/18);
bh = @(V) 1/(exp(-(V+35)/10)+1);

% anonymous functions for channel, leakage, and pump currents
I_K = @(n, V) gK * power(n, 4) * (V - V_K_disp);
I_L = @(V) gL * (V - V_leak);
I_Na = @(V, m, h) gNa * power(m,3) * h * (V - V_Na_disp);
I_P = -I_L(V(1));

dndt = @(an, bn, n, t) an*(1-n) - bn*n;     % rate of change of n in ms^-1
dmdt = @(am, bm, m, t) am*(1-m) - bm*m;     % rate of change of m in ms^-1
dhdt = @(ah, bh, h, t) ah*(1-h) - bh*h;     % rate of change of h in ms^-1

% rage of change of membrane potential
dVdt = @(I_ext, I_K, I_Na, I_L, V, t) (I_ext - I_K - I_Na - I_L - I_P)/C;

Na_finished_open = false;  % to keep track whether Na channel finish or not
Na_open = false;
K_open = false;
for i = 2 : N
    % update t
    t(i) = i * dt;
    
    % determine if there is stimulus based on time
    if (t(i-1) >= stimulus_start) && (t(i-1) <= (stimulus_start + ...
                                                 stimulus_duration))
        I_ext = Iext;
    else
        I_ext = 0;
    end
    
    % compute the coefficients in the equations
    an_0 = an(V(i-1));
    am_0 = am(V(i-1));
    ah_0 = ah(V(i-1));
    bn_0 = bn(V(i-1));
    bm_0 = bm(V(i-1));
    bh_0 = bh(V(i-1));
    
    % using rk4 to determine dn, dm, and dh
    dn1 = dndt(an_0, bn_0, n(i-1), t(i-1)) * dt;
    dn2 = dndt(an_0, bn_0, n(i-1) + 0.5*dn1,...
        t(i-1) + 0.5*dt) * dt;
    dn3 = dndt(an_0, bn_0, n(i-1) + 0.5*dn2,...
        t(i-1) + 0.5*dt) * dt;
    dn4 = dndt(an_0, bn_0, n(i-1) + dn1,...
        t(i-1) + dt) * dt;
    
    dm1 = dmdt(am_0, bm_0, m(i-1), t(i-1)) * dt;
    dm2 = dmdt(am_0, bm_0, m(i-1) + 0.5*dm1,...
        t(i-1) + 0.5*dt) * dt;
    dm3 = dmdt(am_0, bm_0, m(i-1) + 0.5*dm2,...
        t(i-1) + 0.5*dt) * dt;
    dm4 = dmdt(am_0, bm_0, m(i-1) + dm3,...
        t(i-1) + dt) * dt;
    
    dh1 = dhdt(ah_0, bh_0, h(i-1), t(i-1)) * dt;
    dh2 = dhdt(ah_0, bh_0, h(i-1) + 0.5*dh1, t(i-1) + 0.5*dt) * dt;
    dh3 = dhdt(ah_0, bh_0, h(i-1) + 0.5*dh2, t(i-1) + 0.5*dt) * dt;
    dh4 = dhdt(ah_0, bh_0, h(i-1) + dh1, t(i-1) + dt) * dt;
    
    % update n, m, h
    n(i) = n(i-1) +(dn1 + 2*dn2 + 2*dn3 + dn4)/6;
    m(i) = m(i-1) +(dm1 + 2*dm2 + 2*dm3 + dm4)/6;
    h(i) = h(i-1) +(dh1 + 2*dh2 + 2*dh3 + dh4)/6;

    % calculate the existing currents
    if K_o(i-1) >= 0 && Na_i(i-1) >= 0
        I_P = -I_L(V(1));
        dK = dKdt * dt;
        dNa = dNadt * dt;
        K_i(i) = K_i(i-1) - dK;
        K_o(i) = K_o(i-1) + dK;
        Na_i(i) = Na_i(i-1) + dNa;
        Na_o(i) = Na_o(i-1) - dNa;
    else
        I_P = 0;
    end
    
    % using the rule of the model to determine whether there is Na current,
    % K current, or neither of them.
    if V(i-1) >= V_Na_start && Na_finished_open == false
        % when V is larger than the initial condition for Na channel to
        % open and Na channel has not finished opening yet, use the
        % equation    
        Na_open = true;
    else
        Na_open = false;       
    end
    
    if V(i-1) >= V_Na_thres     
        % if V is larger than Na's threshold,
        % then Na channel has finished opening
        Na_finished_open = true;
        K_open = true;
    end
    
    if  Na_finished_open == true && V(i-1) <= V(i-2)
        % if the Na channel has finished opening and V keeps decreasing
        % then K channel is open
        K_open = true;     
    elseif K_open == true && V(i-1) >= V_Na_thres
        % safty measure to make sure that if there is still a postive dv in
        % the current stage, the next step won't close the K channel
        K_open = true;
    else
        K_open = false;
    end
    
    % Na channel open
    if Na_open == true && K_open == false
        if V(i-1) <= -54.4
            I_L_1 = I_L(V(i-1));  
        else
            I_L_1 = 0;
        end
        I_Na_1 = I_Na(V(i-1), m(i-1), h(i-1));
        I_K_1 = 0;
        dv1 = dVdt(I_ext, I_K_1, I_Na_1, I_L_1, V(i-1), t(i-1)) * dt;
        
        if V(i-1) <= -54.4
            I_L_2 = I_L(V(i-1) + 0.5*dv1);  
        else
            I_L_2 = 0;
        end
        I_Na_2 = I_Na(V(i-1) + 0.5*dv1, m(i-1) + 0.5*dm1,...
            h(i-1) + 0.5*dh1);
        I_K_2 = 0;
        dv2 = dVdt(I_ext, I_K_2, I_Na_2, I_L_2, V(i-1) + 0.5*dv1,...
            t(i-1) + 0.5*dt) * dt;

        if V(i-1) <= -54.4
            I_L_3 = I_L(V(i-1) + 0.5*dv2);  
        else
            I_L_3 = 0;
        end
        I_Na_3 = I_Na(V(i-1) + 0.5*dv2, m(i-1) + 0.5*dm2,...
            h(i-1) + 0.5*dh2);
        I_K_3 = 0;
        dv3 = dVdt(I_ext, I_K_3, I_Na_3, I_L_3, V(i-1) + 0.5*dv2,...
            t(i-1) + 0.5*dt) * dt;

        if V(i-1) <= -54.4
            I_L_4 = I_L(V(i-1) + dv3);  
        else
            I_L_4 = 0;
        end
        I_Na_4 = I_Na(V(i-1) + dv3, m(i-1) + dm3, h(i-1)...
            + dh3);
        I_K_4 = 0;
        dv4 = dVdt(I_ext, I_K_4, I_Na_4, I_L_4, V(i-1) + dv3,...
            t(i-1) + dt) * dt;
        
    % K channel open
    elseif K_open == true
        if V(i-1) <= -54.4
            I_L_1 = I_L(V(i-1));  
        else
            I_L_1 = 0;
        end
        I_Na_1 = 0;
        I_K_1 = I_K(n(i-1), V(i-1));
        dv1 = dVdt(I_ext, I_K_1, I_Na_1, I_L_1, V(i-1), t(i-1)) * dt;

        if V(i-1) <= -54.4
            I_L_2 = I_L(V(i-1) + 0.5*dv1);  
        else
            I_L_2 = 0;
        end
        I_Na_2 = 0;
        I_K_2 = I_K(n(i-1) + 0.5*dn1, V(i-1) + 0.5*dv1);
        dv2 = dVdt(I_ext, I_K_2, I_Na_2, I_L_2, V(i-1) + 0.5*dv1,...
            t(i-1) + 0.5*dt) * dt;

        if V(i-1) <= -54.4
            I_L_3 = I_L(V(i-1) + 0.5*dv2);  
        else
            I_L_3 = 0;
        end
        I_Na_3 = 0;
        I_K_3 = I_K(n(i-1) + 0.5*dn2, V(i-1) + 0.5*dv2);
        dv3 = dVdt(I_ext, I_K_3, I_Na_3, I_L_3, V(i-1) + 0.5*dv2,...
            t(i-1) + 0.5*dt) * dt;

        if V(i-1) <= -54.4
            I_L_4 = I_L(V(i-1) + dv3);  
        else
            I_L_4 = 0;
        end
        I_Na_4 = 0;
        I_K_4 = I_K(n(i-1) + dn3, V(i-1) + dv3);
        dv4 = dVdt(I_ext, I_K_4, I_Na_4, I_L_4, V(i-1) + dv3,...
            t(i-1) + dt) * dt;
        
    % before stimulus / after hyperpolarization
    else
        if V(i-1) <= -54.4
            I_L_1 = I_L(V(i-1));  
        else
            I_L_1 = 0;
        end
        I_Na_1 = 0;
        I_K_1 = 0;
        dv1 = dVdt(I_ext, I_K_1, I_Na_1, I_L_1, V(i-1), t(i-1)) * dt;

        if V(i-1) <= -54.4
            I_L_2 = I_L(V(i-1) + 0.5*dv1);  
        else
            I_L_2 = 0;
        end
        I_Na_2 = 0;
        I_K_2 = 0;
        dv2 = dVdt(I_ext, I_K_2, I_Na_2, I_L_2, V(i-1) + 0.5*dv1,...
            t(i-1) + 0.5*dt) * dt;

        if V(i-1) <= -54.4
            I_L_3 = I_L(V(i-1) + 0.5*dv2);  
        else
            I_L_3 = 0;
        end
        I_Na_3 = 0;
        I_K_3 = 0;
        dv3 = dVdt(I_ext, I_K_3, I_Na_3, I_L_3, V(i-1) + 0.5*dv2,...
            t(i-1) + 0.5*dt) * dt;

        if V(i-1) <= -54.4
            I_L_4 = I_L(V(i-1) + dv3);  
        else
            I_L_4 = 0;
        end
        I_Na_4 = 0;
        I_K_4 = 0;
        dv4 = dVdt(I_ext, I_K_4, I_Na_4, I_L_4, V(i-1) + dv3,...
            t(i-1) + dt) * dt;
    end
    
    % update V
    V(i) = V(i-1) +(dv1 + 2*dv2 + 2*dv3 + dv4)/6;    
end

% plot n/m/h v.s. time on the same figure
figure
hold on;
plot(t, n);
plot(t, m);
plot(t, h);
title("Graph for n, m, and h versus time")
legend('n K activation gating','m Na activation gating ','h K inactivation gating ')
xlabel('time(s)')
hold off

% plot potential v.s. time
figure;
plot(t, V);
title("Action potential with membrane potential versus time")
xlabel('time (ms)')
ylabel('memberane potential (mv)')