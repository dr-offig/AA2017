function [SPL, SPL_FF, VR, VT, VRL, VTL, THETA, PRESSURE] = piston(F,a,r,M)
% [SPL, SPL_FARFIELD_APPROX, V_RADIAL, V_TANGENTIAL, V_RADIAL_RMS, V_TANGENTIAL_RMS] = piston(frequency, speaker_radius, distance_from_piston)

  % constants
  c = 340;        % speed of sound in metres per second
  rho = 1.225;    % density of air in kg per cubic metre
  p0 = 20e-6;     % defined as 0dB reference -- 20 microPascals
  v0 = 1;         % indicative piston velocity amplitude. Doesn't really matter for relative polar plots
  N = 40;         % maximum number of terms to calculate in the expansion (recommended value is 40)
  delta = 0.0001; % tolerance for convergence

  % calculated parameters
  omega = F * 2 * pi;
  k = omega / c;
  T = 1 / F;
  
  % calculate coefficients
  ka = k * a;
  A = calculate_recursive_coefficients(ka,N);
  B = calculate_constant_coefficients(ka,N,A);

  % angle and time grid
  Acount = 179;
  THETA = zeros(Acount,1);    % to be filled in with the actual angles of the grid
  PRESSURE = zeros(Acount,1);
  PRESSURE_FF = zeros(Acount,1);       % far field approximation for comparison
  VELOCITY_RADIAL = zeros(Acount,1);
  VELOCITY_TANGENT = zeros(Acount,1);

  for aa = [1:Acount]

    % record the value of theta
    theta = -(pi/2) + (pi/(Acount+1)) + ((aa-1) / (Acount+1)) * pi;
    THETA(aa) = theta;

    % First term is a bit different to the rest (because the derivative of the zeroth hankel function has a different form)
    b = B(1);

    p = b * calculate_pressure_term(theta,k,r,0);
    pressure = p * rho * c; % * exp(-1i * omega * t);

    % vr = ((b * 1i) / (rho * k * c)) * legendreP(0, cos(theta)) * k * (-1) * besselh(1,k*r);
    vr = ((b * 1i) / (rho * k * c)) * calculate_radial_velocity_term(theta,k,r,0);
    velocity_r =  vr; % * exp(-1i * omega * t);

    vt = ((b * 1i) / (rho * k * c)) * calculate_tangent_velocity_term(theta,k,r,0);
    velocity_t = vt; % * exp(-1i * omega * t);

    for n = [1:N-1]
      b = B(n+1);
      p = b * rho * c * calculate_pressure_term(theta,k,r,n);
      pressure = pressure + p; %% * exp(-1i * omega * t);

      vr = ((b * 1i) / (rho * k * c)) * calculate_radial_velocity_term(theta,k,r,n);
      velocity_r = velocity_r + vr; %% * exp(-1i * omega * t);

      vt = ((b * 1i) / (rho * k * c)) * calculate_tangent_velocity_term(theta,k,r,n);
      velocity_t = velocity_t + vt; %% * exp(-1i * omega * t);

      % only keep calculating terms until less than 1 percent improvement
      if ((abs(p) / abs(pressure)) < delta)
        break
      end  

    end
    PRESSURE(aa) = v0 * pressure; % * exp(-1i * omega * t);
    VELOCITY_RADIAL(aa) = v0 * velocity_r; % * exp(-1i * omega * t);
    VELOCITY_TANGENT(aa) = v0 * velocity_t; % * exp(-1i * omega * t);

    % also calculate the FarField approximation for comparison
    if (theta == 0)
      pressure_ff = rho * c * exact_on_axis_pressure(k,r,a);
    else
      pressure_ff = -1i * k * a^2 * rho * c * (besselj(1, ka * sin(theta)) / (ka * sin(theta))) * (exp(1i * k * r) / r);
    end
    PRESSURE_FF(aa) = v0 * pressure_ff; % * exp(-1i * omega * t);

  end  


  %%% Outputs %%%
  SPL = 20 * log10 (abs(PRESSURE)./p0);
  SPL_FF = 20 * log10 (abs(PRESSURE_FF)./p0);
  
  VRL = abs(VELOCITY_RADIAL);
  VTL = abs(VELOCITY_TANGENT);
  VT = harmonic_motion(VELOCITY_TANGENT,zeros(Acount,1),M);
  VR = harmonic_motion(VELOCITY_RADIAL,zeros(Acount,1),M);
  
end

%%%%%%% private functions %%%%%%%
function A = calculate_recursive_coefficients(ka,N)
  
  A = ones(N,1);
  
  A(1) = ka^2 / 2;
  for  n=[1:N-1]
    index = n+1;   %% the recursion formula starts with n = 0, but matlab indexes arrays starting at 1
    A(index) = A(index-1) * (-1) * ka^2 * ( (2*n - 1) / (32 * n^3 - 26 * n + 6));
  end

end


function B = calculate_constant_coefficients(ka,N,A)

  B = ones(size(A));
  
  for n=[0:N-1]
    index = n+1;
    B(index) = A(index) * hypergeom(n+1,[n+2, 2*n + 3/2], -1 * (ka/2)^2);
  end
  
end


function p = exact_on_axis_pressure(k,r,a)
  p = -2i * exp(1i*k*(r + (r^2+a^2)^(1/2))/2) * sin((k*(r^2+a^2)^(1/2) - k*r)/2); 
end


function y = spherical_hankel(n,x)
  y = (pi/(2*x))^(1/2)*besselh(n+1/2,x);
end



function p = calculate_pressure_term(theta,k,r,n)
  p = legendreP(2*n, cos(theta)) * spherical_hankel(2*n,k*r);
end


function vr = calculate_radial_velocity_term(theta,k,r,n)
  %vr = legendreP(2*n, cos(theta)) * (k/2) * (besselh(2*n-1,k*r) - besselh(2*n+1,k*r));
  %vr = legendreP(2*n, cos(theta)) * k * (spherical_hankel(2*n-1,k*r) - (1/(k*r)) * (spherical_hankel(2*n,k*r) + k*r*spherical_hankel(2*n+1,k*r)));
  vr = legendreP(2*n, cos(theta)) * k * ( ((2*n)/(k*r)) * spherical_hankel(2*n,k*r) - spherical_hankel(2*n+1,k*r) );
end


function vt = calculate_tangent_velocity_term(theta,k,r,n)
  if (theta == 0)
    vt = 0;
  else
    vt =  ((2*n+1)/(r * sin(theta))) * (cos(theta) * legendreP(2*n,cos(theta)) - legendreP(2*n+1,cos(theta))) * spherical_hankel(2*n,k*r);
  end
end 
