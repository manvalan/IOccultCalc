# Guida Completa Formato .oop - IOccultCalc Object-Oriented Preset
# =====================================================================
# 
# Il formato .oop è la configurazione standard di IOccultCalc per definire
# tutti i parametri di una sessione di calcolo occultazioni in modo strutturato.
#
# Versione: 2.0 (30 Novembre 2025)
# Autore: Michele Bigi
# Compatibilità: IOccultCalc v2.0+ con AstDyn Integration

# =====================================================================
# STRUTTURA GENERALE
# =====================================================================
#
# Il file .oop è organizzato in SEZIONI, ognuna con un nome seguito da punto:
#
#   sezione_name.
#       .keyword = valore
#       .altro_keyword = altro_valore
#
# REGOLE SINTASSI:
# - I nomi sezione terminano sempre con "."
# - Le keyword iniziano sempre con "."
# - I valori possono essere: numeri, stringhe, booleani, array
# - Commenti: linee che iniziano con "#" o "!"
# - Stringhe: racchiuse in 'apici' o "virgolette"
# - Booleani: .TRUE. / .FALSE. (FORTRAN style)
# - Array: [elemento1, elemento2, elemento3]

# =====================================================================
# SEZIONE GENERAL - Configurazione Globale
# =====================================================================
# STATUS: **OBBLIGATORIA**
# Contiene i parametri principali della sessione di calcolo.

general.
    # PROPAGATORE - Algoritmo di propagazione orbitale
    # OBBLIGATORIO | Valori: 'chebyshev' | 'rkf78' | 'astdyn' | 'hybrid'
    .propagator = 'hybrid'
    
    # STEP SIZE - Passo temporale base [giorni]  
    # OBBLIGATORIO | Range: 0.001 - 10.0 | Default: 0.01
    .step_size_days = 0.01
    
    # TOLERANCE - Tolleranza numerica integratore
    # OPZIONALE | Range: 1e-15 - 1e-6 | Default: 1e-12
    .tolerance = 1e-12
    
    # SEARCH RADIUS - Raggio ricerca [gradi]
    # OBBLIGATORIO | Range: 0.001 - 5.0 | Default: 1.0  
    .search_radius_deg = 1.0
    
    # PRECISION MODE - Livello precisione calcoli
    # OPZIONALE | Valori: 'fast' | 'balanced' | 'precise' | Default: 'balanced'
    .precision_mode = 'precise'
    
    # VERBOSE - Output dettagliato
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .verbose = .FALSE.
    
    # PARALLEL - Calcolo parallelo
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .parallel = .TRUE.
    
    # MAX THREADS - Numero massimo thread
    # OPZIONALE | Range: 1 - 64 | Default: auto-detect
    .max_threads = 8

# =====================================================================
# SEZIONE ORBIT_SOURCE - Origine Elementi Orbitali
# =====================================================================
# STATUS: **OBBLIGATORIA**
# Definisce da dove caricare gli elementi orbitali degli asteroidi.

orbit_source.
    # TYPE - Tipo di sorgente elementi
    # OBBLIGATORIO | Valori: 'mpc' | 'astdys' | 'jpl' | 'file' | 'manual'
    .type = 'astdys'
    
    # FILE - Path file elementi (se type = 'file')
    # CONDIZIONALE | Required se type = 'file' | Formati: .eq1, .oel, .txt
    .file = 'examples/asteroids.txt'
    
    # URL - URL elementi online (se type = 'mpc'/'jpl')
    # CONDIZIONALE | Required se type = 'mpc'/'jpl'
    .url = 'https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT'
    
    # CACHE_DIR - Directory cache elementi
    # OPZIONALE | Default: './cache/elements/'
    .cache_dir = './cache/elements/'
    
    # CACHE_EXPIRE - Scadenza cache [ore]
    # OPZIONALE | Range: 1 - 8760 | Default: 24
    .cache_expire_hours = 24
    
    # UPDATE_AUTO - Aggiornamento automatico
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .update_auto = .TRUE.

# =====================================================================
# SEZIONE ORBIT_FITTING - Controllo Orbit Fitting con Osservazioni
# =====================================================================
# STATUS: **OPZIONALE** 
# Controlla l'uso di osservazioni astrometriche per migliorare gli elementi.

orbit_fitting.
    # MODE - Modalità orbit fitting
    # OBBLIGATORIO | Valori: 'never' | 'auto' | 'always' | Default: 'auto'
    .mode = 'auto'
    
    # OBSERVATIONS_SOURCE - Sorgente osservazioni
    # OPZIONALE | Valori: 'astdys' | 'mpc' | 'file' | Default: 'astdys'
    .observations_source = 'astdys'
    
    # OBSERVATIONS_FILE - File osservazioni locali
    # CONDIZIONALE | Required se observations_source = 'file' | Formato: .rwo
    .observations_file = 'observations/433.rwo'
    
    # MIN_OBSERVATIONS - Minimo osservazioni per fitting
    # OPZIONALE | Range: 3 - 1000 | Default: 5
    .min_observations = 5
    
    # MIN_ARC_DAYS - Minimo arco osservativo [giorni]
    # OPZIONALE | Range: 1 - 36500 | Default: 30
    .min_arc_days = 30
    
    # OUTLIER_THRESHOLD - Soglia outlier detection [sigma]
    # OPZIONALE | Range: 1.0 - 5.0 | Default: 3.0
    .outlier_threshold = 3.0
    
    # MAX_ITERATIONS - Massime iterazioni differential correction
    # OPZIONALE | Range: 1 - 50 | Default: 20
    .max_iterations = 20
    
    # CONVERGENCE_TOLERANCE - Tolleranza convergenza [AU]
    # OPZIONALE | Range: 1e-9 - 1e-3 | Default: 1e-6
    .convergence_tolerance = 1e-6
    
    # FORCE_REFIT - Forza refit ogni predizione
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .force_refit = .FALSE.
    
    # LOG_DETAILS - Log dettagliato fitting
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .log_details = .FALSE.

# =====================================================================
# SEZIONE ASTEROIDS - Selezione Asteroidi Target
# =====================================================================
# STATUS: **OBBLIGATORIA**
# Definisce quali asteroidi includere nei calcoli.

asteroids.
    # SELECTION_MODE - Modalità selezione asteroidi
    # OBBLIGATORIO | Valori: 'all' | 'list' | 'range' | 'magnitude' | 'file'
    .selection_mode = 'list'
    
    # LIST - Lista numeri asteroidi (se selection_mode = 'list')
    # CONDIZIONALE | Array di interi | Required se selection_mode = 'list'
    .list = [433, 1, 2, 4, 10, 15, 16, 18, 29, 324]
    
    # RANGE - Range numeri (se selection_mode = 'range')
    # CONDIZIONALE | Array [min, max] | Required se selection_mode = 'range'
    .range = [1, 1000]
    
    # FILE - File lista asteroidi (se selection_mode = 'file')
    # CONDIZIONALE | Path file | Required se selection_mode = 'file'
    .file = 'examples/priority_asteroids.txt'
    
    # MAG_LIMIT - Limite magnitudine (se selection_mode = 'magnitude')
    # CONDIZIONALE | Range: 5.0 - 25.0 | Required se selection_mode = 'magnitude'
    .mag_limit = 15.0
    
    # DIAMETER_MIN - Diametro minimo [km] (filtro opzionale)
    # OPZIONALE | Range: 0.1 - 1000.0
    .diameter_min = 1.0
    
    # DIAMETER_MAX - Diametro massimo [km] (filtro opzionale)
    # OPZIONALE | Range: 0.1 - 1000.0
    .diameter_max = 100.0
    
    # EXCLUDE - Lista asteroidi da escludere
    # OPZIONALE | Array di interi
    .exclude = [2060, 3200, 4179]
    
    # PRIORITY - Lista asteroidi alta priorità (processati per primi)
    # OPZIONALE | Array di interi
    .priority = [433, 1036, 25143]

# =====================================================================
# SEZIONE TIME - Intervallo Temporale di Ricerca
# =====================================================================
# STATUS: **OBBLIGATORIA**
# Definisce il periodo di tempo per cercare le occultazioni.

time.
    # START_DATE - Data inizio ricerca
    # OBBLIGATORIO | Formato: 'YYYY-MM-DD' o 'YYYY-MM-DD HH:MM:SS'
    .start_date = '2026-01-01'
    
    # END_DATE - Data fine ricerca
    # OBBLIGATORIO | Formato: 'YYYY-MM-DD' o 'YYYY-MM-DD HH:MM:SS'
    .end_date = '2026-12-31'
    
    # TIME_ZONE - Fuso orario output
    # OPZIONALE | Valori: 'UTC' | 'local' | '+HH:MM' | Default: 'UTC'
    .time_zone = 'UTC'
    
    # STEP_SIZE - Passo temporale ricerca [ore]
    # OPZIONALE | Range: 0.1 - 24.0 | Default: 1.0
    .step_size_hours = 1.0
    
    # EPOCH_REFERENCE - Epoca di riferimento elementi
    # OPZIONALE | Formato: 'YYYY-MM-DD' | Default: auto (epoca elementi)
    .epoch_reference = '2026-06-15'

# =====================================================================
# SEZIONE OBSERVER - Configurazione Osservatore
# =====================================================================
# STATUS: **OPZIONALE**
# Definisce la posizione dell'osservatore (default: centro Terra).

observer.
    # LOCATION_TYPE - Tipo posizione osservatore
    # OBBLIGATORIO | Valori: 'geocentric' | 'topocentric' | 'spacecraft'
    .location_type = 'topocentric'
    
    # LONGITUDE - Longitudine [gradi] (se topocentric)
    # CONDIZIONALE | Range: -180.0 - +180.0 | Required se topocentric
    .longitude = 11.2558  # Roma
    
    # LATITUDE - Latitudine [gradi] (se topocentric)
    # CONDIZIONALE | Range: -90.0 - +90.0 | Required se topocentric  
    .latitude = 41.9028   # Roma
    
    # ALTITUDE - Altitudine [metri] (se topocentric)
    # OPZIONALE | Range: -500 - 9000 | Default: 0
    .altitude = 50
    
    # NAME - Nome osservatorio
    # OPZIONALE | Stringa descrittiva
    .name = 'Roma, Italia'
    
    # IAU_CODE - Codice IAU osservatorio
    # OPZIONALE | Codice a 3 caratteri MPC
    .iau_code = '095'

# =====================================================================
# SEZIONE SEARCH - Parametri Ricerca Occultazioni
# =====================================================================
# STATUS: **OBBLIGATORIA**
# Configura i criteri per identificare le occultazioni.

search.
    # MAX_MAGNITUDE - Magnitudine stellare massima
    # OBBLIGATORIO | Range: 5.0 - 20.0 | Default: 12.0
    .max_magnitude = 14.0
    
    # MIN_MAGNITUDE - Magnitudine stellare minima (filtro opzionale)
    # OPZIONALE | Range: 3.0 - 18.0
    .min_magnitude = 8.0
    
    # SEARCH_RADIUS_DEG - Raggio ricerca intorno asteroide [gradi]
    # OBBLIGATORIO | Range: 0.001 - 2.0 | Default: 0.05
    .search_radius_deg = 0.1
    
    # MIN_PROBABILITY - Probabilità minima occultazione [0-1]
    # OPZIONALE | Range: 0.001 - 1.0 | Default: 0.01
    .min_probability = 0.05
    
    # MIN_DURATION - Durata minima [secondi]
    # OPZIONALE | Range: 0.001 - 3600.0 | Default: 0.1
    .min_duration_sec = 0.1
    
    # MAX_DURATION - Durata massima [secondi] (filtro outlier)
    # OPZIONALE | Range: 0.1 - 7200.0 | Default: 600
    .max_duration_sec = 120.0
    
    # MIN_DROP - Drop magnitudine minimo
    # OPZIONALE | Range: 0.01 - 15.0 | Default: 0.1
    .min_drop_mag = 0.3
    
    # SOLAR_ELONGATION_MIN - Elongazione solare minima [gradi]
    # OPZIONALE | Range: 0 - 180 | Default: 15
    .solar_elongation_min = 20.0
    
    # LUNAR_DISTANCE_MIN - Distanza lunare minima [gradi]
    # OPZIONALE | Range: 0 - 180 | Default: 5
    .lunar_distance_min = 10.0
    
    # ALTITUDE_MIN - Altitudine minima [gradi]
    # OPZIONALE | Range: 0 - 90 | Default: 15
    .altitude_min = 20.0

# =====================================================================
# SEZIONE STAR_CATALOGS - Cataloghi Stellari
# =====================================================================
# STATUS: **OPZIONALE**
# Configura i cataloghi stellari da utilizzare.

star_catalogs.
    # PRIMARY - Catalogo stellare principale
    # OBBLIGATORIO | Valori: 'gaia_dr3' | 'gaia_dr2' | 'hipparcos' | 'tycho2' | 'ucac4'
    .primary = 'gaia_dr3'
    
    # SECONDARY - Catalogo di backup
    # OPZIONALE | Valori: 'gaia_dr2' | 'hipparcos' | 'tycho2'
    .secondary = 'gaia_dr2'
    
    # LOCAL_PATH - Path catalogo locale
    # OPZIONALE | Directory con file catalogo
    .local_path = './catalogs/gaia_dr3/'
    
    # DOWNLOAD_AUTO - Download automatico se mancante
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .download_auto = .TRUE.
    
    # CACHE_SIZE_GB - Dimensione massima cache [GB]
    # OPZIONALE | Range: 1 - 100 | Default: 10
    .cache_size_gb = 20
    
    # PROPER_MOTION - Usa moti propri stellari
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .proper_motion = .TRUE.
    
    # PARALLAX_CORRECTION - Correzione parallasse
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .parallax_correction = .TRUE.

# =====================================================================
# SEZIONE ASTDYN - Configurazione AstDyn (se propagator = 'astdyn')
# =====================================================================
# STATUS: **OPZIONALE**
# Parametri specifici per il propagatore AstDyn RKF78.

astdyn.
    # TOLERANCE - Tolleranza integratore RKF78
    # OPZIONALE | Range: 1e-15 - 1e-6 | Default: 1e-12
    .tolerance = 1e-12
    
    # MIN_STEP - Passo minimo integratore [giorni]
    # OPZIONALE | Range: 1e-8 - 1e-3 | Default: 1e-6
    .min_step_days = 1e-6
    
    # MAX_STEP - Passo massimo integratore [giorni]
    # OPZIONALE | Range: 0.1 - 100.0 | Default: 10.0
    .max_step_days = 10.0
    
    # USE_PLANETS - Perturbazioni planetarie
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .use_planets = .TRUE.
    
    # USE_ASTEROIDS - Perturbazioni asteroidali (AST17)
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .use_asteroids = .TRUE.
    
    # USE_RELATIVITY - Correzioni relativistiche
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .use_relativity = .TRUE.
    
    # ASTEROID_PERTURBERS - Lista asteroidi perturbatori
    # OPZIONALE | Array numeri | Default: [1,2,4,10,15,16,18,29,52,511,704,911,1036,1627,2060,3200,4179]
    .asteroid_perturbers = [1, 2, 4, 10, 15, 16]
    
    # COORDINATE_SYSTEM - Sistema coordinate
    # OPZIONALE | Valori: 'ecliptic' | 'equatorial' | Default: 'ecliptic'
    .coordinate_system = 'ecliptic'

# =====================================================================
# SEZIONE CHEBYSHEV - Configurazione Chebyshev (se propagator = 'chebyshev')
# =====================================================================
# STATUS: **OPZIONALE**
# Parametri per i polinomi di Chebyshev per screening veloce.

chebyshev.
    # DEGREE - Grado polinomi Chebyshev
    # OPZIONALE | Range: 4 - 20 | Default: 8
    .degree = 8
    
    # INTERVAL_DAYS - Intervallo validità [giorni]
    # OPZIONALE | Range: 1 - 90 | Default: 30
    .interval_days = 30
    
    # REGENERATE_THRESHOLD - Soglia rigenerazione [giorni]
    # OPZIONALE | Range: 0.1 - 10.0 | Default: 1.0
    .regenerate_threshold_days = 1.0
    
    # ACCURACY_TARGET - Accuratezza target [arcsec]
    # OPZIONALE | Range: 0.001 - 10.0 | Default: 0.1
    .accuracy_target_arcsec = 0.05
    
    # CACHE_POLYNOMIALS - Cache polinomi su disco
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .cache_polynomials = .TRUE.

# =====================================================================
# SEZIONE HYBRID - Configurazione Strategia Ibrida
# =====================================================================
# STATUS: **OPZIONALE** 
# Parametri per strategia ibrida Chebyshev (FASE 1) + RKF78 (FASE 2).

hybrid.
    # PHASE1_ALGORITHM - Algoritmo FASE 1 (screening)
    # OPZIONALE | Valori: 'chebyshev' | 'keplergeocentric' | Default: 'chebyshev'
    .phase1_algorithm = 'chebyshev'
    
    # PHASE2_ALGORITHM - Algoritmo FASE 2 (closest approach)
    # OPZIONALE | Valori: 'rkf78' | 'astdyn' | Default: 'astdyn'
    .phase2_algorithm = 'astdyn'
    
    # PROMOTION_THRESHOLD - Soglia promozione FASE 1→2 [arcsec]
    # OPZIONALE | Range: 10 - 7200 | Default: 60
    .promotion_threshold_arcsec = 30.0
    
    # CONVERGENCE_TOLERANCE - Tolleranza convergenza FASE 2 [arcsec]
    # OPZIONALE | Range: 0.01 - 10.0 | Default: 1.0
    .convergence_tolerance_arcsec = 0.5
    
    # SEARCH_WINDOW - Finestra ricerca closest approach [ore]
    # OPZIONALE | Range: 0.5 - 12.0 | Default: 2.0
    .search_window_hours = 4.0
    
    # CACHE_PHASE1 - Cache risultati FASE 1
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .cache_phase1 = .TRUE.

# =====================================================================
# SEZIONE OUTPUT - Configurazione Output
# =====================================================================
# STATUS: **OBBLIGATORIA**
# Definisce formato e destinazione dei risultati.

output.
    # FORMAT - Formato output principale
    # OBBLIGATORIO | Valori: 'iota' | 'json' | 'csv' | 'xml' | 'txt'
    .format = 'iota'
    
    # DIRECTORY - Directory output
    # OPZIONALE | Default: './output/'
    .directory = './results/'
    
    # FILENAME - Nome file base (senza estensione)
    # OPZIONALE | Default: 'occultations_YYYY-MM-DD'
    .filename = 'occultations_2026'
    
    # INCLUDE_HEADER - Include header descrittivo
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .include_header = .TRUE.
    
    # INCLUDE_STATISTICS - Include statistiche sessione
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .include_statistics = .TRUE.
    
    # INCLUDE_ELEMENTS - Include elementi orbitali usati
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .include_elements = .FALSE.
    
    # INCLUDE_FITTING_INFO - Include info orbit fitting
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .include_fitting_info = .TRUE.
    
    # PRECISION_COORDINATES - Cifre decimali coordinate [arcsec]
    # OPZIONALE | Range: 1 - 6 | Default: 2
    .precision_coordinates = 3
    
    # PRECISION_TIME - Cifre decimali tempi [secondi]  
    # OPZIONALE | Range: 1 - 6 | Default: 1
    .precision_time = 2
    
    # SORT_BY - Criterio ordinamento
    # OPZIONALE | Valori: 'time' | 'probability' | 'magnitude' | 'duration' | Default: 'time'
    .sort_by = 'probability'
    
    # SORT_ORDER - Ordine ordinamento
    # OPZIONALE | Valori: 'asc' | 'desc' | Default: 'asc'
    .sort_order = 'desc'
    
    # FILTER_DUPLICATES - Filtra eventi duplicati
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .filter_duplicates = .TRUE.
    
    # COMPRESSION - Compressione output
    # OPZIONALE | Valori: 'none' | 'gzip' | 'zip' | Default: 'none'
    .compression = 'none'

# =====================================================================
# SEZIONE PERFORMANCE - Ottimizzazione Performance
# =====================================================================
# STATUS: **OPZIONALE**
# Configurazioni per ottimizzare le prestazioni.

performance.
    # PARALLEL_ASTEROIDS - Parallelizzazione per asteroidi
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .parallel_asteroids = .TRUE.
    
    # PARALLEL_STARS - Parallelizzazione per stelle
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .parallel_stars = .TRUE.
    
    # THREAD_POOL_SIZE - Dimensione pool thread
    # OPZIONALE | Range: 1 - 128 | Default: auto-detect
    .thread_pool_size = 16
    
    # MEMORY_LIMIT_GB - Limite memoria [GB]
    # OPZIONALE | Range: 1 - 1024 | Default: 8
    .memory_limit_gb = 16
    
    # CACHE_STRATEGY - Strategia caching
    # OPZIONALE | Valori: 'none' | 'memory' | 'disk' | 'hybrid' | Default: 'hybrid'
    .cache_strategy = 'hybrid'
    
    # PREFETCH_CATALOGS - Pre-caricamento cataloghi
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .prefetch_catalogs = .TRUE.
    
    # BATCH_SIZE - Dimensione batch processamento
    # OPZIONALE | Range: 1 - 10000 | Default: 100
    .batch_size = 500
    
    # PROGRESS_UPDATE - Intervallo aggiornamento progresso [%]
    # OPZIONALE | Range: 1 - 50 | Default: 10
    .progress_update_percent = 5

# =====================================================================
# SEZIONE DEBUG - Configurazione Debug e Log
# =====================================================================
# STATUS: **OPZIONALE**
# Parametri per debugging e logging dettagliato.

debug.
    # LOG_LEVEL - Livello di logging
    # OPZIONALE | Valori: 'error' | 'warning' | 'info' | 'debug' | 'trace' | Default: 'info'
    .log_level = 'debug'
    
    # LOG_FILE - File di log
    # OPZIONALE | Default: stdout
    .log_file = './logs/ioccultcalc.log'
    
    # LOG_PROPAGATION - Log dettagli propagazione
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .log_propagation = .FALSE.
    
    # LOG_FITTING - Log dettagli orbit fitting
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .log_fitting = .TRUE.
    
    # LOG_PERFORMANCE - Log statistiche performance
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .log_performance = .TRUE.
    
    # SAVE_INTERMEDIATE - Salva risultati intermedi
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .save_intermediate = .FALSE.
    
    # DUMP_ELEMENTS - Dump elementi orbitali per debug
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .dump_elements = .FALSE.
    
    # VALIDATE_INPUT - Validazione rigorosa input
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .validate_input = .TRUE.

# =====================================================================
# SEZIONE ADVANCED - Opzioni Avanzate
# =====================================================================
# STATUS: **OPZIONALE**
# Parametri per utenti esperti e casi d'uso specializzati.

advanced.
    # COORDINATE_FRAME - Frame coordinate di riferimento
    # OPZIONALE | Valori: 'icrs' | 'fk5' | 'j2000' | 'b1950' | Default: 'icrs'
    .coordinate_frame = 'icrs'
    
    # TIME_SCALE - Scala temporale
    # OPZIONALE | Valori: 'utc' | 'tdb' | 'tt' | Default: 'tdb'
    .time_scale = 'tdb'
    
    # ABERRATION_CORRECTION - Correzione aberrazione stellare
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .aberration_correction = .TRUE.
    
    # LIGHT_TIME_CORRECTION - Correzione tempo luce
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .TRUE.
    .light_time_correction = .TRUE.
    
    # ATMOSPHERIC_REFRACTION - Rifrazione atmosferica
    # OPZIONALE | Valori: .TRUE. | .FALSE. | Default: .FALSE.
    .atmospheric_refraction = .FALSE.
    
    # NUTATION_MODEL - Modello nutazione
    # OPZIONALE | Valori: 'iau2000' | 'iau1980' | 'none' | Default: 'iau2000'
    .nutation_model = 'iau2000'
    
    # PRECESSION_MODEL - Modello precessione
    # OPZIONALE | Valori: 'iau2000' | 'iau1976' | Default: 'iau2000'
    .precession_model = 'iau2000'
    
    # PLANETARY_EPHEMERIS - Effemeridi planetarie
    # OPZIONALE | Valori: 'de440' | 'de430' | 'de421' | Default: 'de440'
    .planetary_ephemeris = 'de440'
    
    # CUSTOM_PERTURBERS - Perturbatori personalizzati
    # OPZIONALE | Array [number, mass_ratio, a, e, i] per ogni perturbatore
    .custom_perturbers = []
    
    # INTEGRATION_METHOD - Metodo integrazione numerica
    # OPZIONALE | Valori: 'rkf78' | 'dop853' | 'adams' | 'radau' | Default: 'rkf78'
    .integration_method = 'rkf78'
    
    # ERROR_CONTROL - Controllo errori numerici
    # OPZIONALE | Valori: 'strict' | 'balanced' | 'relaxed' | Default: 'balanced'
    .error_control = 'strict'

# =====================================================================
# ESEMPI CONFIGURAZIONI PREDEFINITE
# =====================================================================

# ESEMPIO 1: Configurazione veloce per survey ampi
survey_veloce.
    .propagator = 'hybrid'
    .step_size_days = 0.1
    .search_radius_deg = 0.5
    .precision_mode = 'fast'
    .mode = 'never'  # No orbit fitting
    .max_magnitude = 12.0
    .format = 'csv'

# ESEMPIO 2: Configurazione precisione per occultazioni finali  
precisione_finale.
    .propagator = 'astdyn'
    .step_size_days = 0.001
    .tolerance = 1e-14
    .search_radius_deg = 0.01
    .precision_mode = 'precise'
    .mode = 'always'  # Sempre orbit fitting
    .max_magnitude = 16.0
    .format = 'iota'
    .include_fitting_info = .TRUE.

# ESEMPIO 3: Configurazione bilanciata per uso generale
uso_generale.
    .propagator = 'hybrid'
    .step_size_days = 0.01
    .search_radius_deg = 0.05
    .precision_mode = 'balanced'
    .mode = 'auto'  # Orbit fitting automatico
    .max_magnitude = 14.0
    .format = 'json'
    .parallel = .TRUE.

# =====================================================================
# NOTE FINALI
# =====================================================================
#
# 1. SEZIONI OBBLIGATORIE: general, orbit_source, asteroids, time, search, output
# 2. SEZIONI OPZIONALI: Tutte le altre (parametri hanno valori default)
# 3. VALIDAZIONE: IOccultCalc valida automaticamente tutti i parametri
# 4. BACKWARD COMPATIBILITY: Supporto per formati .oop legacy
# 5. ESTENSIBILITÀ: Nuove sezioni/parametri possono essere aggiunte
#
# Per esempi pratici, vedere:
# - preset_*.oop nella directory principale
# - examples/preset_*.oop per casi d'uso specifici
#
# Documentazione completa: README.md e PRESET_GUIDE.md
#
# =====================================================================
# Fine Guida Completa Formato .oop v2.0
# =====================================================================