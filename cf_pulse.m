%%%%%Climbing Fiber Calcium Pulse%%%%%%%
function cf_pulse(pulse_time)
	global pfpc_model;
	label = num2str(pulse_time*200);
	event_name = ['event' label];
	ca_on = [event_name 'ca_on'];
	ca_off = [event_name 'ca_off'];
	off_time = pulse_time+0.002;
	%Create parameters
	on_param = [event_name 'on'];
	time_on_cond = ['time >= ' on_param];
	on_param = addparameter(pfpc_model,on_param, 'Value', pulse_time, 'ValueUnits', 'second')
	off_param = [event_name 'off'];
	time_off_cond = ['time >= ' off_param];
	off_param = addparameter(pfpc_model,off_param, 'Value', off_time, 'ValueUnits', 'second')
	%create triggers for ca
	ca_on = addevent(pfpc_model, time_on_cond, 'kinflux = ca_pulse_cf')
	ca_off = addevent(pfpc_model, time_off_cond, 'kinflux = 0')
end
