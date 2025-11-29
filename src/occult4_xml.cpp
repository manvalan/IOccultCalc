#include "ioccultcalc/occult4_xml.h"
#include "ioccultcalc/time_utils.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>

namespace ioccultcalc {

// ============================================================================
// Constructor/Destructor
// ============================================================================

Occult4XMLHandler::Occult4XMLHandler() {
    // Inizializza libxml2
    xmlInitParser();
}

Occult4XMLHandler::~Occult4XMLHandler() {
    // Cleanup libxml2
    xmlCleanupParser();
}

// Helper: Converti JD in data calendario
void jdToCalendar(double jd, int& year, int& month, int& day, double& ut) {
    int a = (int)(jd + 0.5);
    int b = a + 1537;
    int c = (int)((b - 122.1) / 365.25);
    int d = (int)(365.25 * c);
    int e = (int)((b - d) / 30.6001);
    
    day = b - d - (int)(30.6001 * e);
    month = e - (e < 14 ? 1 : 13);
    year = c - (month > 2 ? 4716 : 4715);
    
    double frac = jd + 0.5 - (int)(jd + 0.5);
    ut = frac * 24.0;  // UT in decimal hours
}

// ============================================================================
// IMPORT da XML
// ============================================================================

std::vector<Occult4XMLHandler::Occult4Event> 
Occult4XMLHandler::loadFromXML(const std::string& filename) {
    // Leggi file
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open XML file: " + filename);
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string xmlContent = buffer.str();
    file.close();
    
    return parseXML(xmlContent);
}

std::vector<Occult4XMLHandler::Occult4Event>
Occult4XMLHandler::parseXML(const std::string& xmlContent) {
    std::vector<Occult4Event> events;
    
    // Parse XML
    xmlDocPtr doc = xmlReadMemory(xmlContent.c_str(), xmlContent.length(),
                                  "noname.xml", nullptr, 0);
    if (doc == nullptr) {
        throw std::runtime_error("Failed to parse XML content");
    }
    
    // Get root element
    xmlNode* root = xmlDocGetRootElement(doc);
    if (root == nullptr) {
        xmlFreeDoc(doc);
        throw std::runtime_error("Empty XML document");
    }
    
    // Find all <Event> or <Occultation> nodes
    for (xmlNode* node = root->children; node; node = node->next) {
        if (node->type != XML_ELEMENT_NODE) continue;
        
        std::string nodeName = (const char*)node->name;
        if (nodeName == "Event" || nodeName == "Occultation") {
            Occult4Event event = parseEventNode(node);
            events.push_back(event);
        }
    }
    
    xmlFreeDoc(doc);
    return events;
}

Occult4XMLHandler::Occult4Event 
Occult4XMLHandler::parseEventNode(void* nodePtr) {
    xmlNode* node = (xmlNode*)nodePtr;
    Occult4Event event;
    
    // Parse child elements
    for (xmlNode* child = node->children; child; child = child->next) {
        if (child->type != XML_ELEMENT_NODE) continue;
        
        std::string name = (const char*)child->name;
        
        // Asteroid info
        if (name == "AsteroidNumber") {
            event.asteroidNumber = extractTextContent(child);
        } else if (name == "AsteroidName") {
            event.asteroidName = extractTextContent(child);
        } else if (name == "AsteroidDesignation") {
            event.asteroidDesignation = extractTextContent(child);
        } else if (name == "AsteroidMag") {
            event.asteroidMag = extractDoubleContent(child);
        }
        
        // Star info
        else if (name == "StarCatalog") {
            event.starCatalog = extractTextContent(child);
        } else if (name == "StarID") {
            event.starId = extractTextContent(child);
        } else if (name == "StarRA") {
            event.starRA = extractDoubleContent(child);
        } else if (name == "StarDec") {
            event.starDec = extractDoubleContent(child);
        } else if (name == "StarMag") {
            event.starMag = extractDoubleContent(child);
        } else if (name == "StarDistance") {
            event.starDistance = extractDoubleContent(child);
        }
        
        // Event time
        else if (name == "JulianDate" || name == "JD") {
            event.jdEvent = extractDoubleContent(child);
        } else if (name == "DateTime" || name == "EventTime") {
            event.dateTimeUTC = extractTextContent(child);
        }
        
        // Geometry
        else if (name == "CloseApproach" || name == "CA") {
            event.closeApproachDist = extractDoubleContent(child);
        } else if (name == "PositionAngle" || name == "PA") {
            event.posAngle = extractDoubleContent(child);
        } else if (name == "PathWidth") {
            event.pathWidth = extractDoubleContent(child);
        } else if (name == "MaxDuration") {
            event.maxDuration = extractDoubleContent(child);
        } else if (name == "Uncertainty") {
            event.uncertainty = extractDoubleContent(child);
        } else if (name == "Probability") {
            event.probability = extractDoubleContent(child);
        } else if (name == "MagDrop") {
            event.dropMag = extractDoubleContent(child);
        } else if (name == "EventID") {
            event.eventId = extractTextContent(child);
        }
        
        // Path points
        else if (name == "CenterLine" || name == "CentralPath") {
            event.centerLine = parsePathPoints(child);
        } else if (name == "NorthLimit") {
            event.northLimit = parsePathPoints(child);
        } else if (name == "SouthLimit") {
            event.southLimit = parsePathPoints(child);
        }
    }
    
    return event;
}

std::vector<Occult4XMLHandler::Occult4Event::PathPoint>
Occult4XMLHandler::parsePathPoints(void* nodePtr) {
    xmlNode* node = (xmlNode*)nodePtr;
    std::vector<Occult4Event::PathPoint> points;
    
    for (xmlNode* child = node->children; child; child = child->next) {
        if (child->type != XML_ELEMENT_NODE) continue;
        
        std::string name = (const char*)child->name;
        if (name == "Point" || name == "PathPoint") {
            Occult4Event::PathPoint point;
            
            for (xmlNode* attr = child->children; attr; attr = attr->next) {
                if (attr->type != XML_ELEMENT_NODE) continue;
                
                std::string attrName = (const char*)attr->name;
                if (attrName == "Latitude" || attrName == "Lat") {
                    point.latitude = extractDoubleContent(attr);
                } else if (attrName == "Longitude" || attrName == "Lon") {
                    point.longitude = extractDoubleContent(attr);
                } else if (attrName == "JD") {
                    point.jd = extractDoubleContent(attr);
                } else if (attrName == "DateTime") {
                    point.dateTime = extractTextContent(attr);
                } else if (attrName == "Altitude" || attrName == "StarAlt") {
                    point.altitude = extractDoubleContent(attr);
                } else if (attrName == "SunAltitude" || attrName == "SunAlt") {
                    point.sunAltitude = extractDoubleContent(attr);
                }
            }
            
            points.push_back(point);
        }
    }
    
    return points;
}

OccultationEvent Occult4XMLHandler::toIOccultCalcEvent(const Occult4Event& o4) {
    OccultationEvent event;
    
    // Asteroid data
    event.asteroid.name = o4.asteroidName;
    event.asteroid.designation = o4.asteroidDesignation.empty() ? 
                                 o4.asteroidNumber : o4.asteroidDesignation;
    
    // Star data
    event.star.sourceId = o4.starId;
    event.star.pos.ra = o4.starRA * DEG_TO_RAD;
    event.star.pos.dec = o4.starDec * DEG_TO_RAD;
    event.star.phot_g_mean_mag = o4.starMag;
    
    // Event timing
    event.timeCA.jd = o4.jdEvent;
    
    // Geometry
    event.closeApproachDistance = o4.closeApproachDist;
    event.positionAngle = o4.posAngle;
    event.maxDuration = o4.maxDuration;
    event.probability = o4.probability;
    
    // Uncertainty
    event.uncertaintyNorth = o4.uncertainty;
    event.uncertaintySouth = o4.uncertainty;
    
    // Shadow path
    event.shadowPath.clear();
    for (const auto& pt : o4.centerLine) {
        ShadowPathPoint ioPoint;
        ioPoint.location.latitude = pt.latitude * DEG_TO_RAD;
        ioPoint.location.longitude = pt.longitude * DEG_TO_RAD;
        ioPoint.location.altitude = 0.0;
        ioPoint.time.jd = pt.jd;
        ioPoint.duration = 0.0; // non disponibile da XML
        ioPoint.centerlineDistance = 0.0;
        event.shadowPath.push_back(ioPoint);
    }
    
    event.eventId = o4.eventId;
    
    return event;
}

// ============================================================================
// EXPORT a XML
// ============================================================================

bool Occult4XMLHandler::exportToXML(const OccultationEvent& event, 
                                    const std::string& filename) {
    std::string xml = generateXML(event);
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    file << xml;
    file.close();
    return true;
}

bool Occult4XMLHandler::exportMultipleToXML(
    const std::vector<OccultationEvent>& events,
    const std::string& filename) {
    
    std::string xml = generateXML(events);
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    file << xml;
    file.close();
    return true;
}

std::string Occult4XMLHandler::generateXML(const OccultationEvent& event) {
    std::ostringstream xml;
    xml << std::fixed;
    
    xml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml << "<Occelmnt generator=\"IOccultCalc\" version=\"1.0\">\n";
    
    // Metadata
    xml << "  <Metadata>\n";
    xml << "    <Source>" << escapeXML(options_.organizationName) << "</Source>\n";
    
    // Current timestamp
    time_t now = time(nullptr);
    struct tm* timeinfo = gmtime(&now);
    char timestamp[32];
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo);
    std::cerr << "[DEBUG XML SINGLE] now=" << now << " timestamp=" << timestamp << std::endl;
    xml << "    <Created>" << timestamp << "</Created>\n";
    xml << "    <Count>1</Count>\n";
    xml << "  </Metadata>\n\n";
    
    // Convert and generate event
    std::cerr << "[DEBUG XML] Calling toOccult4Event..." << std::endl;
    Occult4Event o4 = toOccult4Event(event);
    std::cerr << "[DEBUG XML] EventID=" << o4.eventId << " PathWidth=" << o4.pathWidth << std::endl;
    xml << generateEventXML(o4);
    
    xml << "</Occelmnt>\n";
    
    return xml.str();
}

std::string Occult4XMLHandler::generateXML(
    const std::vector<OccultationEvent>& events) {
    
    std::ostringstream xml;
    xml << std::fixed;
    
    xml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml << "<Occelmnt generator=\"IOccultCalc\" version=\"1.0\">\n";
    
    // Metadata
    xml << "  <Metadata>\n";
    xml << "    <Source>" << escapeXML(options_.organizationName) << "</Source>\n";
    
    // Current timestamp
    time_t now = time(nullptr);
    struct tm* timeinfo = gmtime(&now);
    char timestamp[32];
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo);
    std::cerr << "[DEBUG XML MULTIPLE] now=" << now << " timestamp=" << timestamp << " events=" << events.size() << std::endl;
    xml << "    <Created>" << timestamp << "</Created>\n";
    xml << "    <Count>" << events.size() << "</Count>\n";
    xml << "  </Metadata>\n\n";
    
    // Generate all events
    for (size_t i = 0; i < events.size(); i++) {
        std::cerr << "[DEBUG XML] Event " << (i+1) << "/" << events.size() << std::endl;
        Occult4Event o4 = toOccult4Event(events[i]);
        std::cerr << "[DEBUG XML]   EventID=" << o4.eventId << " PathWidth=" << o4.pathWidth << std::endl;
        xml << generateEventXML(o4);
    }
    
    xml << "</Occelmnt>\n";
    
    return xml.str();
}

// Genera evento nel formato Occult4 CSV
std::string Occult4XMLHandler::generateOccult4EventXML(const OccultationEvent& event) {
    std::ostringstream xml;
    xml << std::fixed << std::setprecision(8);
    
    // Estrai data e ora dall'evento
    JulianDate jd = event.timeCA;
    int year, month, day;
    double ut;
    jdToCalendar(jd.jd, year, month, day, ut);
    
    xml << "  <Event>\n";
    
    // <Elements> - Besselian elements
    // Format: source, duration, year, month, day, UT, x, y, dX, dY, d2X, d2Y, d3X, d3Y
    xml << "    <Elements>";
    xml << options_.organizationName << " " << formatDateTime(2440587.5 + (std::time(nullptr) / 86400.0));
    xml << "," << std::setprecision(2) << event.maxDuration;
    xml << "," << year << "," << month << "," << day;
    xml << "," << std::setprecision(6) << ut;
    // X, Y coordinates (Earth radii) - approssimazione da coordinate geografiche
    xml << ",0.0,0.0";  // x, y at closest approach
    xml << ",0.0,0.0";  // dX, dY (rate of change)
    xml << ",0.0,0.0";  // d2X, d2Y (2nd order)
    xml << ",0.0,0.0";  // d3X, d3Y (3rd order)
    xml << "</Elements>\n";
    
    // <Earth> - Earth orientation
    // Format: SubstellarLong, SubstellarLat, SubsolarLong, SubsolarLat, JWST
    xml << "    <Earth>";
    if (!event.shadowPath.empty()) {
        xml << std::setprecision(6) << event.shadowPath[0].location.longitude * 57.29577951308232; // rad to deg
        xml << "," << event.shadowPath[0].location.latitude * 57.29577951308232;
    } else {
        xml << "0.0,0.0";
    }
    xml << ",0.0,0.0";  // Subsolar long/lat (would need sun position calculation)
    xml << ",0";  // JWST flag
    xml << "</Earth>\n";
    
    // <Star> - Star information
    // Format: ID, RA(hrs), Dec(deg), Mb, Mv, Mr, dia(mas), double, K2, ApparentRA, ApparentDec, 
    //         MdropV, MdropR, AdjustedForNearby, BrightNearbyCount, TotalNearbyCount
    xml << "    <Star>";
    xml << escapeXML(event.star.sourceId);
    xml << "," << std::setprecision(8) << (event.star.pos.ra / 15.0);  // degrees to hours
    xml << "," << event.star.pos.dec;
    xml << "," << std::setprecision(2) << event.star.phot_g_mean_mag;  // Mb
    xml << "," << event.star.phot_g_mean_mag;  // Mv
    xml << "," << event.star.phot_g_mean_mag;  // Mr
    xml << ",0.0";  // stellar diameter (mas)
    xml << ",0";    // double star code
    xml << ",";     // K2 flag (blank)
    xml << "," << std::setprecision(8) << (event.star.pos.ra / 15.0);  // Apparent RA
    xml << "," << event.star.pos.dec;  // Apparent Dec
    double magDrop = 2.0;  // Stima del calo magnitudine
    xml << "," << std::setprecision(2) << magDrop;  // MdropV
    xml << "," << magDrop;  // MdropR
    xml << ",0,-1,-1";  // AdjustedForNearby, BrightNearbyCount, TotalNearbyCount
    xml << "</Star>\n";
    
    // <Object> - Asteroid/Object information
    // Format: number, name, mag, diameter, distance, #rings, #moons, dRA, dDec, 
    //         Taxonomy, DiamUncertainty, PlanetMoonInShadow, MagV, MagR
    xml << "    <Object>";
    xml << escapeXML(event.asteroid.designation);
    xml << "," << escapeXML(event.asteroid.name);
    xml << "," << std::setprecision(2) << event.asteroid.H;  // H magnitude as proxy
    xml << "," << std::setprecision(1) << event.asteroid.diameter;  // diameter in km
    double distance = event.asteroid.a;  // semimajor axis as approximate distance
    xml << "," << std::setprecision(6) << distance;
    xml << ",0,0";  // #rings, #moons
    xml << ",0.0,0.0";  // dRA, dDec (hourly motion)
    xml << ",";  // Taxonomy (blank if unknown)
    xml << "," << std::setprecision(1) << (event.asteroid.diameter * 0.1);  // 10% uncertainty
    xml << ",0";  // PlanetMoonInShadow
    xml << "," << std::setprecision(2) << event.asteroid.H;  // MagV
    xml << "," << event.asteroid.H;  // MagR
    xml << "</Object>\n";
    
    // <Orbit> - Orbital elements (low precision for plotting)
    // Format: equinox, MA, yr, month, day, peri, node, i, e, a, q, H0, CoeffLogR, G
    xml << "    <Orbit>";
    xml << "2000.0";  // Equinox
    xml << ",0.0";    // Mean Anomaly
    xml << "," << year << "," << month << "," << day;
    xml << ",0.0,0.0,0.0";  // peri, node, i
    xml << ",0.0,0.0,0.0";  // e, a, q
    xml << ",0.0,0.0,0.15"; // H0, CoeffLogR, G
    xml << "</Orbit>\n";
    
    // <Errors> - Prediction uncertainties
    // Format: pathWidthFraction, majorAxis, minorAxis, PA, sigma, basis, RUWE, duplicate, nonGaiaPM, UCAC4PM
    xml << "    <Errors>";
    double uncertKm = (event.uncertaintyNorth + event.uncertaintySouth) / 2.0;
    double pathWidthKm = 200.0;  // Valore tipico, dovrebbe essere calcolato
    xml << std::setprecision(3) << (uncertKm / std::max(1.0, pathWidthKm));
    xml << "," << std::setprecision(2) << uncertKm / 6371.0 * 206265.0;  // km to arcsec
    xml << "," << uncertKm / 6371.0 * 206265.0;  // minor axis
    xml << ",0.0";  // PA
    xml << "," << std::setprecision(2) << uncertKm / 6371.0 * 206265.0;  // 1-sigma
    xml << ",Star+Assumed";  // Error basis
    xml << ",1.0";  // RUWE (default)
    xml << ",0,0,0";  // duplicate, nonGaiaPM, UCAC4PM
    xml << "</Errors>\n";
    
    // <ID> - Unique identifier
    // Format: uniqueID, MJD of prediction
    xml << "    <ID>";
    char dateID[32];
    snprintf(dateID, sizeof(dateID), "%04d%02d%02d", year, month, day);
    xml << dateID << "_" << event.star.sourceId.substr(std::max(0, (int)event.star.sourceId.length() - 6));
    xml << "," << std::setprecision(1) << (jd.jd - 2400000.5);  // MJD of prediction
    xml << "</ID>\n";
    
    xml << "  </Event>\n";
    
    return xml.str();
}

std::string Occult4XMLHandler::generateEventXML(const Occult4Event& event) {
    std::ostringstream xml;
    xml << std::fixed << std::setprecision(6);
    
    xml << "  <Event>\n";
    xml << "    <EventID>" << escapeXML(event.eventId) << "</EventID>\n";
    
    // Asteroid
    if (options_.includeAsteroidData) {
        xml << "    <Asteroid>\n";
        if (!event.asteroidNumber.empty()) {
            xml << "      <Number>" << escapeXML(event.asteroidNumber) << "</Number>\n";
        }
        if (!event.asteroidName.empty()) {
            xml << "      <Name>" << escapeXML(event.asteroidName) << "</Name>\n";
        }
        if (!event.asteroidDesignation.empty()) {
            xml << "      <Designation>" << escapeXML(event.asteroidDesignation) << "</Designation>\n";
        }
        if (event.asteroidMag > 0) {
            xml << "      <Magnitude>" << formatDouble(event.asteroidMag, 2) << "</Magnitude>\n";
        }
        xml << "    </Asteroid>\n";
    }
    
    // Star
    if (options_.includeStarData) {
        xml << "    <Star>\n";
        xml << "      <Catalog>" << escapeXML(event.starCatalog) << "</Catalog>\n";
        xml << "      <ID>" << escapeXML(event.starId) << "</ID>\n";
        xml << "      <RA unit=\"degrees\">" << formatDouble(event.starRA, 8) << "</RA>\n";
        xml << "      <Dec unit=\"degrees\">" << formatDouble(event.starDec, 8) << "</Dec>\n";
        xml << "      <RAFormatted>" << formatRA(event.starRA) << "</RAFormatted>\n";
        xml << "      <DecFormatted>" << formatDec(event.starDec) << "</DecFormatted>\n";
        xml << "      <Magnitude>" << formatDouble(event.starMag, 2) << "</Magnitude>\n";
        if (event.starDistance > 0) {
            xml << "      <Distance unit=\"parsec\">" << formatDouble(event.starDistance, 2) << "</Distance>\n";
        }
        xml << "    </Star>\n";
    }
    
    // Event timing
    xml << "    <Time>\n";
    xml << "      <JulianDate>" << formatDouble(event.jdEvent, 8) << "</JulianDate>\n";
    xml << "      <UTC>" << event.dateTimeUTC << "</UTC>\n";
    xml << "    </Time>\n";
    
    // Geometry
    xml << "    <Geometry>\n";
    xml << "      <CloseApproach unit=\"arcsec\">" << formatDouble(event.closeApproachDist, 4) << "</CloseApproach>\n";
    xml << "      <PositionAngle unit=\"degrees\">" << formatDouble(event.posAngle, 2) << "</PositionAngle>\n";
    xml << "      <PathWidth unit=\"km\">" << formatDouble(event.pathWidth, 2) << "</PathWidth>\n";
    xml << "      <MaxDuration unit=\"seconds\">" << formatDouble(event.maxDuration, 2) << "</MaxDuration>\n";
    if (options_.includeUncertainty && event.uncertainty > 0) {
        xml << "      <Uncertainty unit=\"km\">" << formatDouble(event.uncertainty, 2) << "</Uncertainty>\n";
    }
    xml << "      <Probability>" << formatDouble(event.probability, 4) << "</Probability>\n";
    if (event.dropMag > 0) {
        xml << "      <MagnitudeDrop>" << formatDouble(event.dropMag, 2) << "</MagnitudeDrop>\n";
    }
    xml << "    </Geometry>\n";
    
    // Path points
    if (options_.includePathPoints && !event.centerLine.empty()) {
        xml << "    <CenterLine>\n";
        for (const auto& pt : event.centerLine) {
            xml << "      <Point>\n";
            xml << "        <Latitude>" << formatDouble(pt.latitude, 6) << "</Latitude>\n";
            xml << "        <Longitude>" << formatDouble(pt.longitude, 6) << "</Longitude>\n";
            xml << "        <JD>" << formatDouble(pt.jd, 8) << "</JD>\n";
            xml << "        <DateTime>" << pt.dateTime << "</DateTime>\n";
            xml << "        <StarAltitude>" << formatDouble(pt.altitude, 2) << "</StarAltitude>\n";
            xml << "        <SunAltitude>" << formatDouble(pt.sunAltitude, 2) << "</SunAltitude>\n";
            xml << "      </Point>\n";
        }
        xml << "    </CenterLine>\n";
        
        if (options_.includeUncertainty) {
            if (!event.northLimit.empty()) {
                xml << "    <NorthLimit>\n";
                for (const auto& pt : event.northLimit) {
                    xml << "      <Point>\n";
                    xml << "        <Latitude>" << formatDouble(pt.latitude, 6) << "</Latitude>\n";
                    xml << "        <Longitude>" << formatDouble(pt.longitude, 6) << "</Longitude>\n";
                    xml << "      </Point>\n";
                }
                xml << "    </NorthLimit>\n";
            }
            
            if (!event.southLimit.empty()) {
                xml << "    <SouthLimit>\n";
                for (const auto& pt : event.southLimit) {
                    xml << "      <Point>\n";
                    xml << "        <Latitude>" << formatDouble(pt.latitude, 6) << "</Latitude>\n";
                    xml << "        <Longitude>" << formatDouble(pt.longitude, 6) << "</Longitude>\n";
                    xml << "      </Point>\n";
                }
                xml << "    </SouthLimit>\n";
            }
        }
    }
    
    xml << "  </Event>\n\n";
    
    return xml.str();
}

Occult4XMLHandler::Occult4Event 
Occult4XMLHandler::toOccult4Event(const OccultationEvent& event) {
    Occult4Event o4;
    
    // Asteroid
    o4.asteroidName = event.asteroid.name;
    o4.asteroidDesignation = event.asteroid.designation;
    
    // Extract number if present (e.g., "(433) Eros" -> "433")
    if (!o4.asteroidDesignation.empty() && o4.asteroidDesignation[0] == '(') {
        size_t end = o4.asteroidDesignation.find(')');
        if (end != std::string::npos) {
            o4.asteroidNumber = o4.asteroidDesignation.substr(1, end - 1);
        }
    }
    
    // Star
    o4.starCatalog = options_.useGaiaIds ? "Gaia DR3" : "UCAC4";
    o4.starId = event.star.sourceId;
    o4.starRA = event.star.pos.ra * RAD_TO_DEG;
    o4.starDec = event.star.pos.dec * RAD_TO_DEG;
    o4.starMag = event.star.phot_g_mean_mag;
    
    // Time
    o4.jdEvent = event.timeCA.jd;
    o4.dateTimeUTC = TimeUtils::jdToISO(event.timeCA);
    
    // Geometry
    o4.closeApproachDist = event.closeApproachDistance;
    o4.posAngle = event.positionAngle;
    o4.maxDuration = event.maxDuration;
    o4.probability = event.probability;
    o4.uncertainty = (event.uncertaintyNorth + event.uncertaintySouth) / 2.0;
    
    // Path width from asteroid diameter if available
    // (simplified - in reality depends on distance and geometry)
    o4.pathWidth = event.asteroid.diameter > 0 ? event.asteroid.diameter : 200.0;
    
    // Magnitude drop estimate
    o4.dropMag = 2.0; // Simple estimate based on geometry
    
    // Asteroid magnitude
    o4.asteroidMag = event.asteroid.H;
    
    // Event ID - generate if not present
    if (event.eventId.empty()) {
        // Format: ASTEROID_STAR_YYYYMMDD
        char dateStr[16];
        int year, month, day;
        double ut;
        jdToCalendar(event.timeCA.jd, year, month, day, ut);
        snprintf(dateStr, sizeof(dateStr), "%04d%02d%02d", year, month, day);
        
        o4.eventId = o4.asteroidDesignation + "_" + 
                     o4.starId.substr(std::max(0, (int)o4.starId.length() - 8)) + "_" + 
                     dateStr;
    } else {
        o4.eventId = event.eventId;
    }
    
    // Convert shadow path points
    o4.centerLine = generatePathPoints(event.shadowPath);
    
    // If no centerline points, generate a placeholder at geocentric point
    if (o4.centerLine.empty()) {
        Occult4Event::PathPoint pt;
        // Use star's RA/Dec as approximate geocentric point
        // (This is a simplification - real centerline calculation is complex)
        pt.latitude = event.star.pos.dec * RAD_TO_DEG;
        pt.longitude = event.star.pos.ra * RAD_TO_DEG;
        
        // Clamp to valid geographic ranges
        if (pt.latitude > 90.0) pt.latitude = 90.0;
        if (pt.latitude < -90.0) pt.latitude = -90.0;
        while (pt.longitude > 180.0) pt.longitude -= 360.0;
        while (pt.longitude < -180.0) pt.longitude += 360.0;
        
        pt.jd = event.timeCA.jd;
        pt.dateTime = TimeUtils::jdToISO(event.timeCA);
        pt.altitude = 45.0; // Placeholder
        pt.sunAltitude = -20.0; // Assume night time
        
        o4.centerLine.push_back(pt);
    }
    
    return o4;
}

std::vector<Occult4XMLHandler::Occult4Event::PathPoint>
Occult4XMLHandler::generatePathPoints(const std::vector<ShadowPathPoint>& ioPoints) {
    std::vector<Occult4Event::PathPoint> points;
    
    for (const auto& ioPt : ioPoints) {
        Occult4Event::PathPoint pt;
        pt.latitude = ioPt.location.latitude * RAD_TO_DEG;
        pt.longitude = ioPt.location.longitude * RAD_TO_DEG;
        pt.jd = ioPt.time.jd;
        pt.dateTime = TimeUtils::jdToISO(ioPt.time);
        pt.altitude = 0.0; // non disponibile in ShadowPathPoint
        pt.sunAltitude = 0.0; // non disponibile in ShadowPathPoint
        points.push_back(pt);
    }
    
    return points;
}

// ============================================================================
// CONFIGURATION
// ============================================================================

void Occult4XMLHandler::setOptions(const XMLOptions& options) {
    options_ = options;
}

Occult4XMLHandler::XMLOptions Occult4XMLHandler::getOptions() const {
    return options_;
}

bool Occult4XMLHandler::validateXML(const std::string& filename) {
    try {
        loadFromXML(filename);
        return true;
    } catch (...) {
        return false;
    }
}

std::string Occult4XMLHandler::detectXMLVersion(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return "";
    }
    
    std::string line;
    while (std::getline(file, line)) {
        size_t pos = line.find("version=");
        if (pos != std::string::npos) {
            size_t start = line.find('"', pos);
            size_t end = line.find('"', start + 1);
            if (start != std::string::npos && end != std::string::npos) {
                return line.substr(start + 1, end - start - 1);
            }
        }
    }
    
    return "unknown";
}

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

std::string Occult4XMLHandler::extractTextContent(void* nodePtr) {
    xmlNode* node = (xmlNode*)nodePtr;
    xmlChar* content = xmlNodeGetContent(node);
    if (content == nullptr) {
        return "";
    }
    std::string result = (const char*)content;
    xmlFree(content);
    return result;
}

double Occult4XMLHandler::extractDoubleContent(void* nodePtr, double defaultValue) {
    std::string content = extractTextContent(nodePtr);
    if (content.empty()) {
        return defaultValue;
    }
    try {
        return std::stod(content);
    } catch (...) {
        return defaultValue;
    }
}

std::string Occult4XMLHandler::escapeXML(const std::string& text) {
    std::string result;
    for (char c : text) {
        switch (c) {
            case '<':  result += "&lt;"; break;
            case '>':  result += "&gt;"; break;
            case '&':  result += "&amp;"; break;
            case '\'': result += "&apos;"; break;
            case '"':  result += "&quot;"; break;
            default:   result += c; break;
        }
    }
    return result;
}

std::string Occult4XMLHandler::formatDouble(double value, int precision) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << value;
    return ss.str();
}

std::string Occult4XMLHandler::formatRA(double raDeg) {
    double raHours = raDeg / 15.0;
    int h = (int)raHours;
    double mFrac = (raHours - h) * 60.0;
    int m = (int)mFrac;
    double s = (mFrac - m) * 60.0;
    
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(2) << h << ":"
       << std::setw(2) << m << ":"
       << std::fixed << std::setprecision(3) << std::setw(6) << s;
    return ss.str();
}

std::string Occult4XMLHandler::formatDec(double decDeg) {
    char sign = decDeg >= 0 ? '+' : '-';
    decDeg = std::abs(decDeg);
    int d = (int)decDeg;
    double mFrac = (decDeg - d) * 60.0;
    int m = (int)mFrac;
    double s = (mFrac - m) * 60.0;
    
    std::ostringstream ss;
    ss << sign << std::setfill('0') << std::setw(2) << d << ":"
       << std::setw(2) << m << ":"
       << std::fixed << std::setprecision(2) << std::setw(5) << s;
    return ss.str();
}

std::string Occult4XMLHandler::formatDateTime(double jd) {
    JulianDate julianDate;
    julianDate.jd = jd;
    return TimeUtils::jdToISO(julianDate);
}

} // namespace ioccultcalc
