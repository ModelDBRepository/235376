%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                           %
%            			Model of Cerebellar PF-PC 						    %
%	       			Long-Term Depression & Potentiation                     %
%                                                                           %
%                     developed by Andrew R. Gallimore                      %  
%          						December 2017      						    %	
%							for Matlab Simbiology	                        %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This code runs a standard LTD/LTP (prompted to choose) protocol and plots the data

LTD_LTP_Model

prompt = 'LTP? (LTD = 1, LTP = 0) ';
p_type = input(prompt);

if (p_type == 1)
	%%%%LTD Events%%%%%
	%eventObj_NO = addevent(pfpc_model, 'time >= 10000', 'kinflux_no = 0.01')
	%%%%%%PF PULSES%%%%%%
	%%%%%%PF PULSES%%%%%%
	ca_pulse_cf = addparameter(pfpc_model,'ca_pulse_cf', 'Value', 1500, 'ValueUnits', 'micromole/second', 'ConstantValue', false);%2000
	num_pulse = 100;
	pulse_times = [1:1:num_pulse]+10000;
	cf_pulse_times = [1:1:num_pulse]+10000.1;
	for i = 1:num_pulse
		pulse_time = pulse_times(i);
		pf_pulse_ltp(pulse_time)
		disp(i)
	end
	for i = 1:num_pulse
		pulse_time = cf_pulse_times(i);
		cf_pulse(pulse_time)
		disp(i)
	end

%eventObj_NOblock = addevent(pfpc_model, 'time >= 10600', 'kcat_nos = 0')

elseif (p_type == 0)
	%%%%%%300 PF PULSES ALONE%%%%%%
	num_pulse = 300;
	pulse_times = [1:1:num_pulse]+10000;
	for i = 1:num_pulse
		pulse_time = pulse_times(i);
		pf_pulse_ltp(pulse_time)
		disp(i)
	end
else
	disp('Please choose either LTD(1) or LTP(0)')
	clear all
	return
end

simDataObj = sbiosimulate(pfpc_model);

[tSim, AMPARPSD_value] = selectbyname(simDataObj,'ampar_psd');
[tSim, Ca_value] = selectbyname(simDataObj,'Ca');
[tSim, ERK_value] = selectbyname(simDataObj,'ERK_act');
[tSim, CAMKII_value] = selectbyname(simDataObj,'CaMKII_Auton');
[tSim, PKC_value] = selectbyname(simDataObj,'PKC_active');

if (p_type == 1)
	figure(1)
	subplot(2,2,1)       % add first plot in 2 x 2 grid
	plot(tSim, Ca_value, 'r')
	xlim([9900 10200])
	ylim([0 5])
	title('Ca')


	subplot(2,2,2)       % add second plot in 2 x 2 grid
	plot(tSim, AMPARPSD_value, 'b')
	title('AMPARPSD')
	xlim([9000 11800])
	ylim([0 1.5])
	title('AMPAR PSD')

	subplot(2,2,3)       % add third plot in 2 x 2 grid
	plot(tSim, ERK_value, 'b')
	xlim([8000 16000])
	ylim([0 1])
	title('Active ERK')

	subplot(2,2,4)       % add fourth plot in 2 x 2 grid
	plot(tSim, CAMKII_value, 'b')
	xlim([8000 16000])
	ylim([0 100])
	title('Autonomous CaMKII')
	hold on

else
	figure(1)
	subplot(1,2,1)       % add first plot in 2 x 2 grid
	plot(tSim, Ca_value, 'r')
	xlim([9900 10400])
	ylim([0 0.4])
	title('Ca')


	subplot(1,2,2)       % add second plot in 2 x 2 grid
	plot(tSim, AMPARPSD_value, 'r')
	title('AMPARPSD')
	xlim([9000 15000])
	ylim([0 1.5])
	title('AMPAR PSD')
	hold on
end
