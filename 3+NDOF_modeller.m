close all
clear all

%  USER MUST SET TWO VALUES TO DETERMINE PLOT OPTIONS

% Set value for plottype as follows:
%       1 -> use 'ezplot' to plot a graph on normal axis (in rad/s)
%       2 -> produce a semilog plot (in Hz)
%       3 -> produce a graph on linear axis using 'plot' (in Hz)
plottype=3;

% Set value to be positive to produce plots of the modeshapes
modeshape_visualisation = 0;

target_frequencies_floor1 = [];
target_frequencies_floor2 = [];
target_frequencies_floor3 = [];

damping_rate = 0.3;

m = 1.83;                                   % mass of one floor
L = 0.2;                                    % length
b = 0.08;                                   % width
E = 210E9;                                  % Young's Modulus
d = 0.001;                                  % thickness
I = b*d*d*d/12;                             % second moment of area
k = (24*E*I)/(L*L*L);                       % static stiffness for each floor

N = 3 + length(target_frequencies_floor1) + length(target_frequencies_floor2) + length(target_frequencies_floor3);         % number of degrees of freedom = 3 + # mass dampers

total_mass_m = m/5;
m_i = total_mass_m/(N-3);

m_mass_dampers = ones(1,N-3)*m_i;           % set the masses of individual mass dampers
k_mass_dampers_floor1 = ((2*pi*target_frequencies_floor1).^2)*m_i*(1+damping_rate*1j);
k_mass_dampers_floor2 = ((2*pi*target_frequencies_floor2).^2)*m_i*(1+damping_rate*1j);
k_mass_dampers_floor3 = ((2*pi*target_frequencies_floor3).^2)*m_i*(1+damping_rate*1j);

one_all_zeroes = [1;zeros(N-1,1)];          % create a 1xN vector [1;0;0;0;0;.....] 
all_m = [m m m,m_mass_dampers];             % create a Nx1 row vector with [m1 m2 m3 mT1 mT2 ... mTN]

M = diag(all_m);                            % create the mass matrix as a diagonalisation of the above all_m

K = zeros(N,N);                             % create the stiffness matrix, there is probably a more efficient method

K(1,1)=  2*k+sum(k_mass_dampers_floor1); K(2,1)= -1*k;                            K(3,1)=   0 ;
K(1,2)= -1*k;                            K(2,2)=  2*k+sum(k_mass_dampers_floor2); K(3,2)= -1*k;
K(1,3)=   0 ;                            K(2,3)= -1*k;                            K(3,3)=  1*k+sum(k_mass_dampers_floor3);

first_row_K       = [0 0 0,-1*k_mass_dampers_floor1, zeros(1,length(k_mass_dampers_floor2)+length(k_mass_dampers_floor3))];
first_row_K_full  = [first_row_K; zeros(N-1,N)];

second_row_K      = [0 0 0, zeros(1,length(k_mass_dampers_floor1)), -1*k_mass_dampers_floor2, zeros(1,length(k_mass_dampers_floor3))];
second_row_K_full = [zeros(1,N); second_row_K; zeros(N-2,N)];

third_row_K       = [0 0 0, zeros(1,length(k_mass_dampers_floor1)+length(k_mass_dampers_floor2)), -1*k_mass_dampers_floor3];
third_row_K_full  = [zeros(2,N);third_row_K;zeros(N-3,N)];

first_column_K_full  = transpose(first_row_K_full);
second_column_K_full  = transpose(second_row_K_full);
third_column_K_full  = transpose(third_row_K_full);


diagonal1=[0 0 0, k_mass_dampers_floor1, k_mass_dampers_floor2, k_mass_dampers_floor3];
diagonal_of_K = diag(diagonal1);

temp_K = K + first_row_K_full + second_row_K_full + third_row_K_full + first_column_K_full + second_column_K_full + third_column_K_full + diagonal_of_K;
K = temp_K;                                                            

[V,D] = eig(K,M);
syms w;

for imode=1:N
  freqs(imode) = sqrt(D(imode,imode));
end

%  Print natural frequencies and mode vectors in command window
hertz = freqs/(2*pi)
modeshapes = V

B = K - ((w*w)*M); 
% harmonic solution for unit force at floor 1
disp = (inv(B))*one_all_zeroes;

%start of ezplot section
if (plottype == 1)
  hold on
  
  ifloor=1;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line'),'Color','k')
  ifloor=2;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line','Color','b'),'Color','g')
  ifloor=3;
  ezplot(disp(ifloor), [0, 130]);
  set(findobj('Type','line','Color','b'),'Color','r')
  
  set(findobj('Type','line','Color','k'),'Color','b')
  
  set(findobj('Type','line'),'LineStyle','-')
end

% Calculate frequency response functions
all_disp = [];
for w = 1:130;
  B = K - ((w*w)*M); 
  % harmonic solution for unit force at floor 1
  disp = (inv(B))*one_all_zeroes;
  all_disp = [all_disp disp];
end
w = 1:130;

% Log plot
if (plottype == 2)
  semilogy((w./(2*pi)),abs(all_disp(1:3,:)),'-');
  xlabel('Frequency (Hz)')
  ylabel('log Displacement (m)')
  legend('Floor 1','Floor 2', 'Floor 3');

% Linear plot
elseif (plottype == 3)
  plot((w./(2*pi)),(all_disp(1:3,:)),'-');
  xlabel('Frequency (Hz)');
  ylabel('Displacement (m)');
  legend('Floor 1','Floor 2', 'Floor 3');
  %xlim([2 5]);
  %ylim([-0.01 0.01]);
end

% Plot modeshapes

if (modeshape_visualisation > 0 )
  V = [zeros(1,N); V];
  V_ = V + 0.25;
  V = V - 0.25;
  for imode=1:3
    figure
    axis([-5 5 0 3.5])
    title1 = ['Mode ' int2str(imode)];
    title(title1)
    hold on
    plot((V(:,imode)),([0 1 2 3]))
    plot([0 0 0 0],[0 1 2 3],'k')
    plot((V_(:,imode)),([0 1 2 3]))
    for jmode=1:3
      plot([V(jmode+1,imode) V_(jmode+1,imode)],[jmode jmode])
    end
  end
end
