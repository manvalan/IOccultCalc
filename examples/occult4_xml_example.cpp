#include <iostream>
#include <iomanip>
#include <ioccultcalc/occultation_predictor.h>
#include <ioccultcalc/occult4_xml.h>
#include <ioccultcalc/time_utils.h>

using namespace ioccultcalc;

void printEvent(const Occult4XMLHandler::Occult4Event& event) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Event ID: " << event.eventId << "\n";
    std::cout << std::string(70, '-') << "\n";
    
    // Asteroid info
    std::cout << "\nASTEROID:\n";
    if (!event.asteroidNumber.empty()) {
        std::cout << "  Number: " << event.asteroidNumber << "\n";
    }
    if (!event.asteroidName.empty()) {
        std::cout << "  Name: " << event.asteroidName << "\n";
    }
    std::cout << "  Designation: " << event.asteroidDesignation << "\n";
    if (event.asteroidMag > 0) {
        std::cout << "  Magnitude: " << std::fixed << std::setprecision(1) 
                  << event.asteroidMag << "\n";
    }
    
    // Star info
    std::cout << "\nSTAR:\n";
    std::cout << "  Catalog: " << event.starCatalog << "\n";
    std::cout << "  ID: " << event.starId << "\n";
    std::cout << "  RA: " << std::fixed << std::setprecision(6) 
              << event.starRA << "° (J2000)\n";
    std::cout << "  Dec: " << event.starDec << "° (J2000)\n";
    std::cout << "  Magnitude: " << std::fixed << std::setprecision(2) 
              << event.starMag << "\n";
    
    // Event time
    std::cout << "\nEVENT TIME:\n";
    std::cout << "  JD: " << std::fixed << std::setprecision(6) 
              << event.jdEvent << "\n";
    std::cout << "  UTC: " << event.dateTimeUTC << "\n";
    
    // Geometry
    std::cout << "\nGEOMETRY:\n";
    std::cout << "  Close Approach: " << std::fixed << std::setprecision(3) 
              << event.closeApproachDist << " arcsec\n";
    std::cout << "  Position Angle: " << std::fixed << std::setprecision(1) 
              << event.posAngle << "°\n";
    std::cout << "  Path Width: " << std::fixed << std::setprecision(1) 
              << event.pathWidth << " km\n";
    std::cout << "  Max Duration: " << std::fixed << std::setprecision(2) 
              << event.maxDuration << " seconds\n";
    std::cout << "  Uncertainty: ±" << std::fixed << std::setprecision(1) 
              << event.uncertainty << " km\n";
    std::cout << "  Probability: " << std::fixed << std::setprecision(1) 
              << (event.probability * 100) << "%\n";
    
    // Path points
    if (!event.centerLine.empty()) {
        std::cout << "\nCENTRAL PATH: " << event.centerLine.size() << " points\n";
        std::cout << "  First point: " 
                  << std::fixed << std::setprecision(4)
                  << event.centerLine.front().latitude << "°N, "
                  << event.centerLine.front().longitude << "°E\n";
        std::cout << "  Last point: "
                  << event.centerLine.back().latitude << "°N, "
                  << event.centerLine.back().longitude << "°E\n";
    }
    
    std::cout << std::string(70, '=') << "\n";
}

int main(int argc, char* argv[]) {
    std::cout << "IOccultCalc - Occult4 XML Handler Example\n";
    std::cout << "==========================================\n\n";
    
    try {
        // Mode 1: Export predictions to Occult4 XML
        if (argc > 1 && std::string(argv[1]) == "export") {
            std::cout << "Mode: EXPORT to Occult4 XML\n\n";
            
            std::string asteroidDesig = argc > 2 ? argv[2] : "433";
            std::string startDate = argc > 3 ? argv[3] : "2026-01-01";
            std::string endDate = argc > 4 ? argv[4] : "2026-06-30";
            std::string outputFile = argc > 5 ? argv[5] : "occultations.xml";
            
            std::cout << "Searching occultations for asteroid " << asteroidDesig << "\n";
            std::cout << "Period: " << startDate << " to " << endDate << "\n\n";
            
            // Create predictor
            OccultationPredictor predictor;
            predictor.loadAsteroidFromAstDyS(asteroidDesig);
            predictor.setAsteroidDiameter(16.8); // km for (433) Eros
            
            // Find occultations
            JulianDate jdStart = TimeUtils::isoToJD(startDate);
            JulianDate jdEnd = TimeUtils::isoToJD(endDate);
            
            auto events = predictor.findOccultations(
                jdStart, jdEnd,
                13.0,  // mag limit
                0.05,  // search radius
                0.005  // min probability
            );
            
            std::cout << "Found " << events.size() << " occultation events\n\n";
            
            if (events.empty()) {
                std::cout << "No events to export\n";
                return 0;
            }
            
            // Export to Occult4 XML
            Occult4XMLHandler xmlHandler;
            
            // Configure export options
            Occult4XMLHandler::XMLOptions options;
            options.includeUncertainty = true;
            options.includePathPoints = true;
            options.pathPointsResolution = 100;  // Point every 100 km
            options.includeStarData = true;
            options.includeAsteroidData = true;
            options.useGaiaIds = true;
            options.organizationName = "IOccultCalc Project";
            options.observerName = "Example User";
            xmlHandler.setOptions(options);
            
            bool success = xmlHandler.exportMultipleToXML(events, outputFile);
            
            if (success) {
                std::cout << "✓ Successfully exported to: " << outputFile << "\n";
                std::cout << "\nFile can be opened in Occult4 for visualization\n";
            } else {
                std::cerr << "✗ Failed to export XML\n";
                return 1;
            }
        }
        
        // Mode 2: Import from Occult4 XML
        else if (argc > 1 && std::string(argv[1]) == "import") {
            std::cout << "Mode: IMPORT from Occult4 XML\n\n";
            
            if (argc < 3) {
                std::cerr << "Usage: " << argv[0] << " import <xml_file>\n";
                return 1;
            }
            
            std::string xmlFile = argv[2];
            
            std::cout << "Loading XML file: " << xmlFile << "\n\n";
            
            // Create XML handler
            Occult4XMLHandler xmlHandler;
            
            // Detect version
            std::string version = xmlHandler.detectXMLVersion(xmlFile);
            std::cout << "Detected XML version: " << version << "\n\n";
            
            // Load events
            auto events = xmlHandler.loadFromXML(xmlFile);
            
            std::cout << "Loaded " << events.size() << " events from XML\n";
            
            // Print each event
            for (const auto& event : events) {
                printEvent(event);
                
                // Ask if user wants to continue (only for first few)
                if (&event != &events.back() && &event - &events[0] < 5) {
                    std::cout << "\nPress Enter to continue to next event...";
                    std::cin.get();
                }
            }
            
            // Convert to IOccultCalc format
            std::cout << "\n\nConverting to IOccultCalc format...\n";
            std::vector<OccultationEvent> ioEvents;
            for (const auto& o4Event : events) {
                ioEvents.push_back(xmlHandler.toIOccultCalcEvent(o4Event));
            }
            
            std::cout << "✓ Converted " << ioEvents.size() << " events\n";
            std::cout << "\nEvents can now be processed with IOccultCalc tools\n";
        }
        
        // Mode 3: Validate XML file
        else if (argc > 1 && std::string(argv[1]) == "validate") {
            std::cout << "Mode: VALIDATE Occult4 XML\n\n";
            
            if (argc < 3) {
                std::cerr << "Usage: " << argv[0] << " validate <xml_file>\n";
                return 1;
            }
            
            std::string xmlFile = argv[2];
            
            std::cout << "Validating XML file: " << xmlFile << "\n\n";
            
            Occult4XMLHandler xmlHandler;
            
            // Check if file is valid
            bool isValid = xmlHandler.validateXML(xmlFile);
            
            if (isValid) {
                std::cout << "✓ XML file is valid\n";
                
                // Show version and basic info
                std::string version = xmlHandler.detectXMLVersion(xmlFile);
                std::cout << "  Version: " << version << "\n";
                
                auto events = xmlHandler.loadFromXML(xmlFile);
                std::cout << "  Events: " << events.size() << "\n";
                
                return 0;
            } else {
                std::cerr << "✗ XML file is NOT valid\n";
                return 1;
            }
        }
        
        // Default: Show usage
        else {
            std::cout << "Occult4 XML Handler - Import/Export occultation predictions\n\n";
            std::cout << "USAGE:\n\n";
            
            std::cout << "1. Export IOccultCalc predictions to Occult4 XML:\n";
            std::cout << "   " << argv[0] << " export <asteroid> <start_date> <end_date> <output_file>\n";
            std::cout << "   Example: " << argv[0] << " export 433 2026-01-01 2026-06-30 eros.xml\n\n";
            
            std::cout << "2. Import Occult4 XML file:\n";
            std::cout << "   " << argv[0] << " import <xml_file>\n";
            std::cout << "   Example: " << argv[0] << " import predictions.xml\n\n";
            
            std::cout << "3. Validate Occult4 XML file:\n";
            std::cout << "   " << argv[0] << " validate <xml_file>\n";
            std::cout << "   Example: " << argv[0] << " validate occultations.xml\n\n";
            
            std::cout << "FEATURES:\n";
            std::cout << "  • Full Occult4 XML format compatibility\n";
            std::cout << "  • Import predictions from Occult4/Steve Preston\n";
            std::cout << "  • Export IOccultCalc predictions for Occult4\n";
            std::cout << "  • Bidirectional conversion IOccultCalc ↔ Occult4\n";
            std::cout << "  • Gaia DR3 and UCAC4 star catalog support\n";
            std::cout << "  • Shadow path with uncertainty bands\n";
            std::cout << "  • Complete event metadata (time, geometry, probability)\n\n";
        }
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
