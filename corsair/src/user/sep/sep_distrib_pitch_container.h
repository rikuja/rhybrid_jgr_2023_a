/** This file is part of Corsair simulation.
 *
 *  Copyright 2011-2013 Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SEP_DISTRIB_PITCH_CONTAINER_H
#define SEP_DISTRIB_PITCH_CONTAINER_H

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>

#include <map>

namespace sep {

   typedef bool (*finalizePitchDistrib)();
   typedef Real (*getPitchFunction)();
   typedef bool (*initializePitchDistrib)(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName);
   
   struct PitchDistribution {
      finalizePitchDistrib finalize;
      getPitchFunction getPitch;
      initializePitchDistrib initialize;
      
      PitchDistribution();
      PitchDistribution(finalizePitchDistrib finalize,getPitchFunction getPitch,initializePitchDistrib initialize);
   };
   
   class DistribPitchContainer {
    public:
      DistribPitchContainer();
      ~DistribPitchContainer();
      
      bool getDistribution(const std::string& name,finalizePitchDistrib& finalize,
			   getPitchFunction& getPitch,initializePitchDistrib& initialize);
      bool registerDistribution(const std::string& name,finalizePitchDistrib finalize,
				getPitchFunction getPitch,initializePitchDistrib initialize);
    private:
      std::map<std::string,PitchDistribution> pitchDistribs;
	
   };

} // namespace sep

#endif
