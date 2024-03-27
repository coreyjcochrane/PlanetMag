% Evaluate default magnetic field model with no magnetopause for spacecraft flybys

% Part of the PlanetMag framework for evaluation and study of planetary magnetic fields.
% Created by Corey J. Cochrane and Marshall J. Styczinski
% Maintained by Marshall J. Styczinski
% Contact: corey.j.cochrane@jpl.nasa.gov
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WriteTimeSeries("Galileo", 0, -1, 'Io', -2)
WriteTimeSeries("Galileo", 0, -1, 'Europa', -2)
WriteTimeSeries("Galileo", 0, -1, 'Ganymede', -2)
WriteTimeSeries("Galileo", 0, -1, 'Europa', -1, 'modelDescrip', 1)

WriteTimeSeries("Juno", 0, -1, 'Europa')
WriteTimeSeries("Juno", 0, -1, 'Ganymede')

WriteTimeSeries("Cassini", 0, -1, 'Mimas')
