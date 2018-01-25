%%%%%Parallel Fiber Glutamate and Calcium Pulse%%%%%%%
function pf_pulse_ltp(pulse_time)
	global pfpc_model;
	label = num2str(pulse_time*200);
	event_name = ['event' label];
	glu_on = [event_name 'glu_on'];
	glu_off = [event_name 'glu_off'];
	ca_on = [event_name 'ca_on'];
	ca_off = [event_name 'ca_off'];
	off_time = pulse_time+0.001;
	%Create parameters
	on_param = [event_name 'on'];
	time_on_cond = ['time >= ' on_param];
	on_param = addparameter(pfpc_model,on_param, 'Value', pulse_time, 'ValueUnits', 'second')
	off_param = [event_name 'off'];
	time_off_cond = ['time >= ' off_param];
	off_param = addparameter(pfpc_model,off_param, 'Value', off_time, 'ValueUnits', 'second')
	%create triggers for glu and ca
	glu_on = addevent(pfpc_model, time_on_cond, 'kinflux_glu = glu_pulse')
	glu_off = addevent(pfpc_model, time_off_cond, 'kinflux_glu = 0')
	ca_on = addevent(pfpc_model, time_on_cond, 'kinflux_pf = ca_pulse')
	ca_off = addevent(pfpc_model, time_off_cond, 'kinflux_pf = 0')
end
