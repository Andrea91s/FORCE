function  MOrbit_printout(particle, EQLST, sfb)
% function used to printout the results from MORBIT to:
% sfb = 0   screen
% sfb = 1   ASCII file 
% Example; MOrbit_printout(particle, EQLST, sfb = 0)

% Define the output stream
if (sfb == 0)
    fid = 1;
    more off
elseif (sfb == 1)
       filename = ['ORBIT_' num2str(EQLST.pulseNumber) '_at_' num2str(EQLST.selectedtime, '%1.3f') '_s_' ...
                   num2str(particle.initialconditions.energy_eV/1000, '%2.0f') '_keV_' ...
                   num2str(180*particle.initialconditions.pitch_angle/pi, '%3.0f') '_deg_' ...
                   num2str(180*particle.initialconditions.phi/pi, '%3.0f') '_deg_' ... 
                   'R' num2str(particle.orbit.position.R(1), '%1.2f') ...
                   '_Z' num2str(particle.orbit.position.z(1), '%1.2f') '.txt'];
       fid = fopen (filename, "w");
endif

%

fprintf(fid, '----------------------------------------------------------------------------\n')
fprintf(fid, 'ORBIT FOLLOWING CODE\n')
fprintf(fid, '----------------------------------------------------------------------------\n')
fprintf(fid, 'Particle:\n')
fprintf(fid, '\tMass:\t%f (amu)\n', particle.mass/physical_constant ('atomic mass constant'));    
fprintf(fid, '\tCharge:\t%f (e)\n\n', particle.charge/physical_constant ('elementary charge'));

fprintf(fid, 'Initial properties:\n')
fprintf(fid, '\tPosition (x,y,z) =\t(%f, %f, %f) (m)\n', particle.initialconditions.x0, particle.initialconditions.y0, particle.initialconditions.z0);
fprintf(fid, '\tVelocity (vx,vy,vz) =\t(%1.4e, %1.4e, %1.4e) (m/s)\n', particle.initialconditions.vx0, particle.initialconditions.vy0, particle.initialconditions.vz0);

fprintf(fid, '\tVelocity =\t\t%1.4e (m/s)\n', particle.initialconditions.v);
fprintf(fid, '\tVel. parallel =\t\t%1.4e (m/s)\n', particle.initialconditions.v_parl);
fprintf(fid, '\tVel. perperndicular =\t%1.4e (m/s)\n', particle.initialconditions.v_perp);
fprintf(fid, '\tEnergy =\t\t%f (eV)\n', particle.initialconditions.energy_eV);
fprintf(fid, '\tMagnetic moment =\t%f (eV/T)\n', energy_unit (particle.orbit.magnetic_moment.mu(1), conversion = 0));
fprintf(fid, '\tPitch angle =\t\t%f (rad)\n', particle.initialconditions.pitch_angle);
fprintf(fid, '\tLambda =\t\t%f\n', particle.initialconditions.lambda);
fprintf(fid, '\tGyro-angle =\t\t%f (rad)\n', particle.initialconditions.phi);
fprintf(fid, '\tLarmor radius =\t\t%f (m)\n', particle.initialconditions.larmor_radius);
fprintf(fid, '\tGyro-frequency =\t%g (Hz), \t%g (rad/s)\n', particle.initialconditions.frequency, particle.initialconditions.omega);

fprintf(fid, 'Fields at initial position:\n')
fprintf(fid, '\tB-field (Bx,By,Bz) =\t(%f, %f, %f) (T)\n', particle.orbit.magnetic_field.Bx(1), particle.orbit.magnetic_field.By(1), particle.orbit.magnetic_field.Bz(1));
fprintf(fid, '\tB-field modulus =\t%f (T)\n', particle.orbit.magnetic_field.B(1));
fprintf(fid, '\tE-field (Ex,Ey,Ez) =\t(%f, %f, %f) (V/m)\n', particle.orbit.electric_field.Ex(1), particle.orbit.electric_field.Ey(1), particle.orbit.electric_field.Ez(1));

fprintf(fid, '\nFinal properties:\n')
fprintf(fid, '\tPosition (x,y,z) =\t(%f, %f, %f) (m)\n', particle.orbit.position.x(end), particle.orbit.position.y(end), particle.orbit.position.z(end));
fprintf(fid, '\tVelocity (vx,vy,vz) =\t(%1.4e, %1.4e, %1.4e) (m/s)\n', particle.orbit.velocity.vx(end), particle.orbit.velocity.vy(end), particle.orbit.velocity.vz(end));
fprintf(fid, '\tVelocity =\t\t%1.4e (m/s)\n', particle.orbit.velocity.v(end));
fprintf(fid, '\tVel. parallel =\t\t%1.4e (m/s)\n', particle.orbit.velocity.v_pll(end));
fprintf(fid, '\tVel. perperndicular =\t%1.4e (m/s)\n', particle.orbit.velocity.v_perp(end));
fprintf(fid, '\tEnergy =\t\t%f (eV)\n', energy_unit (particle.orbit.energy.kinetic(end), conversion = 0) );
fprintf(fid, '\tMagnetic moment =\t%f (eV/T)\n', energy_unit (particle.orbit.magnetic_moment.mu(end), conversion = 0));
fprintf(fid, '\tPitch angle =\t\t%f (rad)\n', particle.orbit.pitch_angle(end));
fprintf(fid, '\tLambda =\t\t%f\n', particle.orbit.lambda(end));
fprintf(fid, '\tGyro-angle =\t\t%f (rad)\n', particle.orbit.position.phi(end));
fprintf(fid, '\tLarmor radius =\t\t%f (m)\n', particle.orbit.larmor_radius(end));
fprintf(fid, '\tGyro-frequency =\t%g (Hz) \t%g (rad/s)\n', particle.orbit.angular_frequency(end), particle.orbit.angular_velocity(end));

fprintf(fid, '\nOrbit properties:\n')
fprintf(fid, '\t<LAMBDA>: \t\t%f\n', mean(particle.orbit.LAMBDA));
fprintf(fid, '\t<P/PSIa>: \t\t%f\n', mean(particle.orbit.PtorPSIa));
%fprintf(fid, '\tConfinement: \t\t%s\n', particle.orbit.confinement);
%fprintf(fid, '\tPassing or trapped: \t%s\n', particle.orbit.turning_points.status);
%fprintf(fid, '\tReflection: \t\t%s\n', particle.orbit.vtor_reflection.status);
%fprintf(fid, '\tGyro-ave reflection: \t%s\n', particle.guidingcenter.vtor_reflection.status);
%fprintf(fid, '\tEncircling: \t\t%s\n\n', particle.orbit.encircling.status);

fprintf(fid, 'Fields at final position:\n')
fprintf(fid, '\tB-field (Bx,By,Bz) =\t(%f, %f, %f) (T)\n', particle.orbit.magnetic_field.Bx(end), particle.orbit.magnetic_field.By(end), particle.orbit.magnetic_field.Bz(end));
fprintf(fid, '\tE-field (Ex,Ey,Ez) =\t(%f, %f, %f) (V/m)\n', particle.orbit.electric_field.Ex(end), particle.orbit.electric_field.Ey(end), particle.orbit.electric_field.Ez(end));


fprintf(fid, '----------------------------------------------------------------------------\n')
fprintf(fid, 'GUIDING CENTRE from ORBIT FOLLOWING CODE\n')
fprintf(fid, '----------------------------------------------------------------------------\n')
fprintf(fid, 'Initial properties:\n')
fprintf(fid, '\tPosition (x,y,z) =\t(%f, %f, %f) (m)\n', particle.orbit.guidingcenter_x(1), particle.orbit.guidingcenter_y(1), particle.orbit.guidingcenter_z(1));
fprintf(fid, '\tVelocity (vx,vy,vz) =\t(%1.4e, %1.4e, %1.4e) (m/s)\n', particle.orbit.guidingcenter_VX(1), particle.orbit.guidingcenter_VY(1), particle.orbit.guidingcenter_VZ(1));
fprintf(fid, '\tVelocity =\t\t%1.4e (m/s)\n', particle.orbit.guidingcenter_V(1));
%fprintf(fid, '\tVel. parallel =\t\t%1.4e (m/s)\n', particle.orbit.guidingcenter_v_parl(1));
%fprintf(fid, '\tVel. perperndicular =\t%1.4e (m/s)\n', particle.orbit.guidingcenter_v_perp(1));
%fprintf(fid, '\tEnergy =\t\t%f (eV)\n', energy_unit (particle.orbit.guidingcenter_kinetic_energy(1), conversion = 0));
%fprintf(fid, '\tMagnetic moment =\t%f (eV/T)\n', energy_unit (particle.orbit.guidingcenter_mu(1), conversion = 0));
%fprintf(fid, '\tPitch angle =\t\t%f (rad)\n', particle.orbit.guidingcenter_pitch_angle(1));
%fprintf(fid, '\tLambda =\t\t%f\n', particle.orbit.guidingcenter_lambda(1));

fprintf(fid, 'Fields at initial position:\n')
fprintf(fid, '\tB-field (Bx,By,Bz) =\t(%f, %f, %f) (T)\n', particle.orbit.guidingcenter_Bx(1), particle.orbit.guidingcenter_By(1), particle.orbit.guidingcenter_Bz(1));
fprintf(fid, '\tB-field modulus =\t%f (T)\n', particle.orbit.guidingcenter_B(1));

fprintf(fid, 'Final properties:\n')
fprintf(fid, '\tPosition (x,y,z) =\t(%f, %f, %f) (m)\n', particle.orbit.guidingcenter_x(end), particle.orbit.guidingcenter_y(end), particle.orbit.guidingcenter_z(end));
fprintf(fid, '\tVelocity (vx,vy,vz) =\t(%1.4e, %1.4e, %1.4e) (m/s)\n', particle.orbit.guidingcenter_vx(end), particle.orbit.guidingcenter_vy(end), particle.orbit.guidingcenter_vz(end));
fprintf(fid, '\tVelocity =\t\t%1.4e (m/s)\n', particle.orbit.guidingcenter_V(end));
%fprintf(fid, '\tVelocity =\t\t%1.4e (m/s)\n', particle.guidingcenter.gc.velocity.v(end));
%fprintf(fid, '\tVel. parallel =\t\t%1.4e (m/s)\n', particle.orbit.guidingcenter_v_parl(end));
%fprintf(fid, '\tVel. perperndicular =\t%1.4e (m/s)\n', particle.orbit.guidingcenter_v_perp(end));
%fprintf(fid, '\tEnergy =\t\t%f (eV)\n', energy_unit (particle.guidingcenter.gc.energy.kinetic(end), conversion = 0) );
%fprintf(fid, '\tMagnetic moment =\t%f (eV/T)\n', energy_unit (particle.orbit.magnetic_moment.mu(1), conversion = 0));
%fprintf(fid, '\tPitch angle =\t\t%f (rad)\n', particle.orbit.pitch_angle(end));
%fprintf(fid, '\tLambda =\t\t%f\n', particle.initialconditions.lambda);


if (particle.guidingcenter.gc.solve == 1)
fprintf(fid, '----------------------------------------------------------------------------\n')
fprintf(fid, 'GUIDING CENTRE CODE\n')
fprintf(fid, '----------------------------------------------------------------------------\n')
fprintf(fid, 'Initial properties:\n')
fprintf(fid, '\tPosition (x,y,z) =\t(%f, %f, %f) (m)\n', particle.guidingcenter.gc.position.x0, particle.guidingcenter.gc.position.y0, particle.guidingcenter.gc.position.z0);
%fprintf(fid, '\tVelocity (vx,vy,vz) =\t(%1.4e, %1.4e, %1.4e) (m/s)\n', particle.initialconditions.vx0, particle.initialconditions.vy0, particle.initialconditions.vz0);
%fprintf(fid, '\tVelocity =\t\t%1.4e (m/s)\n', particle.guidingcenter.gc.position.v0);
fprintf(fid, '\tVel. parallel =\t\t%1.4e (m/s)\n', particle.guidingcenter.gc.position.v_par0);
%fprintf(fid, '\tVel. perperndicular =\t%1.4e (m/s)\n', particle.guidingcenter.gc.position.v_per0);
%fprintf(fid, '\tEnergy =\t\t%f (eV)\n', particle.guidingcenter.gc.position.E0);
%fprintf(fid, '\tMagnetic moment =\t%f (eV/T)\n', energy_unit (particle.guidingcenter.gc.position.mu0, conversion = 0));
%fprintf(fid, '\tPitch angle =\t\t%f (rad)\n', particle.initialconditions.pitch_angle);
%fprintf(fid, '\tLambda =\t\t%f\n', particle.initialconditions.lambda);

fprintf(fid, 'Fields at initial position:\n')
fprintf(fid, '\tB-field (Bx,By,Bz) =\t(%f, %f, %f) (T)\n', particle.guidingcenter.gc.Bx0, particle.guidingcenter.gc.By0,particle.guidingcenter.gc.Bz0);
fprintf(fid, '\tB-field modulus =\t%f (T)\n', particle.guidingcenter.gc.B0);

fprintf(fid, 'Final properties:\n')
fprintf(fid, '\tPosition (x,y,z) =\t(%f, %f, %f) (m)\n', particle.guidingcenter.gc.position.x(end), particle.guidingcenter.gc.position.y(end), particle.guidingcenter.gc.position.z(end));
fprintf(fid, '\tVelocity (vx,vy,vz) =\t(%1.4e, %1.4e, %1.4e) (m/s)\n', particle.orbit.velocity.vx(end), particle.orbit.velocity.vy(end), particle.orbit.velocity.vz(end));
%fprintf(fid, '\tVelocity =\t\t%1.4e (m/s)\n', particle.guidingcenter.gc.velocity.v(end));
%fprintf(fid, '\tVel. parallel =\t\t%1.4e (m/s)\n', particle.guidingcenter.gc.velocity.v_par(end));
%fprintf(fid, '\tVel. perperndicular =\t%1.4e (m/s)\n', particle.guidingcenter.gc.velocity.v_per(end));
%fprintf(fid, '\tEnergy =\t\t%f (eV)\n', energy_unit (particle.guidingcenter.gc.energy.kinetic(end), conversion = 0) );
%fprintf(fid, '\tMagnetic moment =\t%f (eV/T)\n', energy_unit (particle.orbit.magnetic_moment.mu(1), conversion = 0));
%fprintf(fid, '\tPitch angle =\t\t%f (rad)\n', particle.orbit.pitch_angle(end));
%fprintf(fid, '\tLambda =\t\t%f\n', particle.initialconditions.lambda);
endif


fprintf(fid, '----------------------------------------------------------------------------\n');
fprintf(fid, 'SIMULATION PARAMETERS:\n');
fprintf(fid, '----------------------------------------------------------------------------\n');
fprintf(fid, '\tstart time:\t\t%e s\n', particle.times.start_time);
fprintf(fid, '\tend time:\t\t%e s\n', particle.times.end_time);
fprintf(fid, '\tsteps:\t\t\t%d\n', particle.times.steps);
fprintf(fid, '\tstep size:\t\t%e s\n', particle.times.dt);
fprintf(fid, '\tmax distance:\t\t%e m\n', particle.orbit.max_distance);

fprintf(fid, 'ODE solver: %s\n', particle.orbit.solver.method);
fprintf(fid, '\tabsolute tolerance:\t%e\n', particle.orbit.solver.abs_tol);
fprintf(fid, '\trelative tolerance:\t%e\n', particle.orbit.solver.rel_tol);
fprintf(fid, '\tintegration method:\t%s\n', particle.orbit.solver.int_method);
fprintf(fid, '\tinitial step size:\t%f\n', particle.orbit.solver.init_step_size);
fprintf(fid, '\tmaximum order:\t\t%f\n', particle.orbit.solver.max_order);
fprintf(fid, '\tmaximum step size:\t%f\n', particle.orbit.solver.max_step_size);
fprintf(fid, '\tminimum step size:\t%f\n', particle.orbit.solver.min_step_size);
fprintf(fid, '\tstep limit :\t\t%f\n', particle.orbit.solver.step_limit);
fprintf(fid, '\tComputation time FO:\t%f (s)\n', particle.orbit.cpu_time);
%if (particle.guidingcenter.gc.solve == 1)
%fprintf(fid, '\tComputation time GC:\t%f (s)\n', particle.guidingcenter.cpu_time);
%endif
fprintf(fid, 'LSODE performance:\n')
fprintf(fid, '\tISTATE: %d\n', particle.orbit.ISTATE);
fprintf(fid, '\tMessage: %s\n', particle.orbit.msg);
fprintf(fid, '\tVelocity change:\t%e (relative)\n', 1 - particle.orbit.velocity.v(end)/particle.orbit.velocity.v(1));
fprintf(fid, '\tMagnetic moment change: %e (relative)\n', 1 - particle.orbit.magnetic_moment.mu(end)/particle.orbit.magnetic_moment.mu(1));
fprintf(fid, '----------------------------------------------------------------------------\n');
     

% Close the stream
if (sfb == 0)
   more off
elseif (sfb == 1)
    fclose (fid);
endif
clear fid
