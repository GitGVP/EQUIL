function pressure = variational_pressure_output( ...
        Pi,Pr,PR,PB,PBB,PrB,PRB,PrR,PRR,Pperp)
%VARIATIONAL_PRESSURE_OUTPUT Named second-order EOS derivative interface.

    pressure.Pi = Pi;
    pressure.Pr = Pr;
    pressure.PR = PR;
    pressure.PB = PB;
    pressure.PBB = PBB;
    pressure.PrB = PrB;
    pressure.PRB = PRB;
    pressure.PrR = PrR;
    pressure.PRR = PRR;
    pressure.Pperp = Pperp;
end
