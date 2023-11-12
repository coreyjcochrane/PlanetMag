function [MagModel, CsheetModel, MPmodel, magModelDescrip, fEnd] = GetModelOpts(parentName, opt, MPopt)
    switch parentName
        case 'Earth'
            % opt 1
            IGRF13 = 'MagFldEarthIGRF13';
            IGRF13sheet = 'None';
            
            if opt == 0; opt = 1; end % Set default to IGRF13
            switch opt
                case 1
                    MagModel = IGRF13;
                    CsheetModel = IGRF13sheet;
                    magModelDescrip = 'IGRF13';
                    fEnd = 'IGRF13';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end

            MPmodel = 'None'; % No magnetopause model implemented
            MPend = 'noMP';

        case 'Jupiter'
            % opt 1
            VIP4 = 'MagFldJupiterVIP4';
            C1981sheet = 'Connerney1981';
            % opt 2
            O6 = 'MagFldJupiterGSFCO6';
            K1997sheet = 'Khurana1997';
            % opt 3
            KhuranaJup = 'VIP4 with O6 orientation';
            K2005sheet = 'KS2005';
            % opt 4
            JRM09 = 'MagFldJupiterJRM09';
            C2020sheet = 'Connerney2020';
            % opt 5 is Vance et al. 2021 combo, JRM09 + C1981
            % opt 6 is Seufert et al. 2011 combo, VIP4 + K1997
            % opt 7 (includes C2020 current sheet)
            JRM33 = 'MagFldJupiterJRM33';
            
            % MPopt 1: Alexeev and Belenkaya (2005) magnetopause model
            AB2005 = 'AB2005'; 
            % MPopt 2: Engle (1992) no-tilt magnetopause field model
            E1992a90 = 'Engle1992alpha90';
            % MPopt 3: Engle (1992) sunward-tilt magnetopause field model
            E1992a0 = 'Engle1992alpha0'; 
            % MPopt 4: Engle (1992) anti-sunward-tilt magnetopause field model
            E1992a180 = 'Engle1992alpha180';
            % MPopt 5: Engle (1992) magnetopause field model with Bode (1994) time-dependent coefficients
            B1994 = 'Bode1994';

            if opt == 0; opt = 7; end % Set default to JRM33 + C2020
            switch opt
                case 1
                    MagModel = VIP4;
                    CsheetModel = C1981sheet;
                    magModelDescrip = 'VIP4 + C1981';
                    fEnd = 'VIP4C1981';
                case 2
                    MagModel = O6;
                    CsheetModel = K1997sheet;
                    magModelDescrip = 'O6 + K1997';
                    fEnd = 'O6K1997';
                case 3
                    MagModel = KhuranaJup;
                    CsheetModel = K2005sheet;
                    magModelDescrip = 'KS2005';
                    fEnd = 'KS2005';
                case 4
                    MagModel = JRM09;
                    CsheetModel = C2020sheet;
                    magModelDescrip = 'JRM09 + C2020';
                    fEnd = 'JRM09C2020';
                case 5
                    MagModel = JRM09;
                    CsheetModel = C1981sheet;
                    magModelDescrip = 'JRM09 + C1981';
                    fEnd = 'JRM09C1981';
                case 6
                    MagModel = VIP4;
                    CsheetModel = K1997sheet;
                    magModelDescrip = 'VIP4 + K1997';
                    fEnd = 'VIP4K1997';
                case 7
                    MagModel = JRM33;
                    CsheetModel = C2020sheet;
                    magModelDescrip = 'JRM33 + C2020';
                    fEnd = 'JRM33C2020';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            if MPopt == 0; MPopt = 2; end % Set default to Engle (1992) no-tilt model
            switch MPopt
                case 1
                    MPmodel = AB2005;
                    MPend = 'AB2005';
                case 2
                    MPmodel = E1992a90;
                    MPend = 'E1992a90';
                case 3
                    MPmodel = E1992a0;
                    MPend = 'E1992a0';
                case 4
                    MPmodel = E1992a180;
                    MPend = 'E1992a180';
                case 5
                    MPmodel = B1994;
                    MPend = 'B1994';
                otherwise
                    MPmodel = 'None';
                    MPend = 'noMP';
            end
            
            
        case 'Saturn'
            % opt 1
            B2010 = 'MagFldSaturnBurton2010';
            B2010sheet = 'SphericalHarmonic';
            % opt 2
            Cassini11 = 'MagFldSaturnCassini11';
            Cassini11sheet = 'Cassini11';

            if opt == 0; opt = 2; end % Set default to Cassini 11
            switch opt
                case 1
                    MagModel = B2010;
                    CsheetModel = B2010sheet;
                    magModelDescrip = 'Burton et al. 2010';
                    fEnd = 'B2010';
                case 2
                    MagModel = Cassini11;
                    CsheetModel = Cassini11sheet;
                    magModelDescrip = 'Cassini 11 field + sheet';
                    fEnd = 'Cassini11';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            MPmodel = 'None'; % No magnetopause model implemented
            MPend = 'noMP';
            
            
        case 'Uranus'
            % opt 1
            Q3 = 'MagFldUranusQ3';
            Q3sheet = 'SphericalHarmonic';
            % opt 2
            AH5 = 'MagFldUranusAH5';
            AH5sheet = 'None';
            
            % MPopt 1: Dipole-mirror magnetopause model used in Herbert
            % (2009) to derive AH5 model (based on SM1996)
            Q3mp = 'Q3mp';
            % MPopt 2: Box harmonic magnetopause model from Arridge and
            % Eggington (2022)
            AE2022 = 'AE2022';

            if opt == 0; opt = 2; end % Set default to AH5
            switch opt
                case 1
                    MagModel = Q3;
                    CsheetModel = 'None'; % Only use uniform external field if noMP
                    magModelDescrip = 'Q3';
                    fEnd = 'Q3';
                case 2
                    MagModel = AH5;
                    CsheetModel = AH5sheet;
                    magModelDescrip = 'AH5';
                    fEnd = 'AH5';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            if MPopt == 0; MPopt = 2; end % Set default to AE2022 model
            switch MPopt
                case 1
                    MPmodel = Q3mp;
                    MPend = 'Q3mp';
                case 2
                    MPmodel = AE2022;
                    MPend = 'AE2022';
                otherwise
                    MPmodel = 'None';
                    MPend = 'noMP';
                    if opt == 1
                        CsheetModel = Q3sheet;
                    end
            end
            
            
        case 'Neptune'
            % opt 1
            O8 = 'MagFldNeptuneO8';
            O8sheet = 'None';
            
            % MPopt 1: Dipole-mirror magnetopause model described
            % in Schulz and McNab (1996)
            SM1996 = 'SM1996';
            
            if opt == 0; opt = 1; end % Set default to O8
            switch opt
                case 1
                    MagModel = O8;
                    CsheetModel = O8sheet;
                    magModelDescrip = 'O8';
                    fEnd = 'O8';
                otherwise
                    warning(['Magnetic field model option ' num2str(opt) ' not recognized. Defaulting to "None".'])
                    MagModel = 'None';
                    CsheetModel = 'None';
                    magModelDescrip = 'None';
                    fEnd = 'None';
            end
            
            if MPopt == 0; MPopt = 1; end % Set default to SM1996 model
            switch MPopt
                case 1
                    MPmodel = SM1996;
                    MPend = 'SM1996';
                otherwise
                    MPmodel = 'None';
                    MPend = 'noMP';
            end
            
    end
    
    if ~strcmp(magModelDescrip, 'None')
        magModelDescrip = [magModelDescrip ' + ' MPend];
        fEnd = [fEnd MPend];
    end
end