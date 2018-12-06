# Soslow
Codebase for Actimetry Core in the Soslow study.
  
  
## Guidelines
  
* Wear Detection
  * use 15 second level files to detect wear status
  * use vector magnitude
* Valid Day
  * During 0500 - 2359, at least 600 minutes of wearing
* Sleep Detection
  * aggregate to 60 second level and then run sleep algorithm
  * map the results back to 15 second level data
  * use axis1
  * use settings for youth
* Physical Activity
  * use 15 second level files to mark activity levels
  * use vector magnitude
  * though cutpoints are for 60 second level observations, divide it by 4 to fit 15 second files
  
## File Structure
```
| /scripts/
   | process_data.R
   | graphs.R
   | tables.R
   | statistical_analysis.R
```
  
## R Packages
* PhysicalActivity 0.2-2
* PhysActBedRest 1.0
* RSQLite 2.1.1
* ggplot2 3.1.0

