clc; clear;

highFidelityPropagator(0, 0, 0, 10, 'atmosphere', 'none', 'maxdeg', 1)

% Function to propagate an orbit using a high-fidelity propagator
% 
% Author: Griffin Jourda 4/01/2023
% 
% Inputs 
%	r0		:	initial position (m) 
%	v0		:	initial velocity (m/s) 
%	jd0		:	epoch Julian date
%	tspan	:	time span to return data for (s)
% Name-value inputs 
%	'atmosphere'	:	atmosphere model, options: 'harrispriester', 'none'
%	'options'		:	ode solver options
%	'maxdeg'		:	maximum degree of gravity model 
%	'eop'			:	earth orientation parameter handling, options:
%						'iau' (default) calculates or interpolates P, N,
%						Pi, dut1, and dtai at every time step;
%						'average' calculates P, N, Pi dut1 and dtai
%						halfway through the simulation and keeps them
%						constant;
%						'none' calculates ECEF vectors via a rotation of
%						GMST only 
% Outputs
%	t	:	simulation time span (s) 
%	y	:	spacecraft states (m, m/s)

function [] = highFidelityPropagator(r0, v0, jd0, tspan, varargin)
	% Process inputs
	if mod(length(varargin), 2) ~= 0
		error('Expected even number of name-value pairs')
	end

	options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12, 'InitialStep', 60, 'MaxStep', 2700);
	model = 'none';
	maxdeg = 0;

	for i = 1:2:length(varargin)
		parameter = varargin{i};
		value = varargin{i+1};

		switch lower(parameter)
			% Atmosphere model options
			case 'atmosphere'
				enforceclass(parameter, value, 'char');
				if ~ismember(lower(value), {'harrispriester', 'none'})
					error("Unrecognized atmosphere option '%s'. Available options are: 'harrispriester', 'none'")
				end

			% ODE Solver options
			case 'options'
				enforceclass(parameter, value, 'struct');
				options = value;
		
			% Gravity model maximum degree
			case 'maxdeg' 
				if ~isnumeric(value) || value < 0 || mod(value, 1) ~= 0
					error('Expected positive integer value for maximum gravity degree')
				elseif value > 0
					[C, S] = loadSHData('EGM2008.txt', maxdeg+1);
				end
			
			% Earth orientation parameter options
			case 'eop'
				enforceclass(parameter, value, 'char')
				if ~ismember(lower(value), {'iau', 'average', 'none'})
					error("Unrecognized EOP option '%s'. Available options are: 'iau', 'average', 'none'")
				end

				switch lower(value)
					case 'iau'
						
					case 'average'
						
					case 'none'
				end
		end
	end

	disp(options) 
	disp(model) 
	disp(maxdeg)
end

function [] = enforceclass(parameter, value, expectedclass)
	if ~strcmpi(class(value), expectedclass)
		error('Expected class %s for parameter %s', expectedclass, parameter);
	end
end