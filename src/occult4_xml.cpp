#include "ioccultcalc/occult4_xml.h"
#include "ioccultcalc/time_utils.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <fstream>
#include <sstream>
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
    Occult4Event o4 = toOccult4Event(event);
    
    std::ostringstream xml;
    xml << std::fixed << std::setprecision(6);
    
    xml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml << "<OccultationPrediction version=\"1.0\" generator=\"IOccultCalc\">\n";
    xml << "  <Metadata>\n";
    xml << "    <Organization>" << escapeXML(options_.organizationName) << "</Organization>\n";
    if (!options_.observerName.empty()) {
        xml << "    <Observer>" << escapeXML(options_.observerName) << "</Observer>\n";
    }
    // Current date/time
    JulianDate now;
    now.jd = 2400000.5 + (std::time(nullptr) / 86400.0);
    xml << "    <GenerationDate>" << formatDateTime(now.jd) << "</GenerationDate>\n";
    xml << "  </Metadata>\n\n";
    
    xml << generateEventXML(o4);
    
    xml << "</OccultationPrediction>\n";
    
    return xml.str();
}

std::string Occult4XMLHandler::generateXML(
    const std::vector<OccultationEvent>& events) {
    
    std::ostringstream xml;
    xml << std::fixed << std::setprecision(6);
    
    xml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml << "<OccultationPrediction version=\"1.0\" generator=\"IOccultCalc\">\n";
    xml << "  <Metadata>\n";
    xml << "    <Organization>" << escapeXML(options_.organizationName) << "</Organization>\n";
    if (!options_.observerName.empty()) {
        xml << "    <Observer>" << escapeXML(options_.observerName) << "</Observer>\n";
    }
    // Current date/time
    JulianDate now;
    now.jd = 2400000.5 + (std::time(nullptr) / 86400.0);
    xml << "    <GenerationDate>" << formatDateTime(now.jd) << "</GenerationDate>\n";
    xml << "    <EventCount>" << events.size() << "</EventCount>\n";
    xml << "  </Metadata>\n\n";
    
    for (const auto& event : events) {
        Occult4Event o4 = toOccult4Event(event);
        xml << generateEventXML(o4);
    }
    
    xml << "</OccultationPrediction>\n";
    
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
    o4.pathWidth = event.asteroid.diameter;
    
    // Event ID
    o4.eventId = event.eventId;
    
    // Convert shadow path points
    o4.centerLine = generatePathPoints(event.shadowPath);
    
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
