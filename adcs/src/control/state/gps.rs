//! A submodule for doing all things gps, holding the 
//! values used for the sidereal time calculations critical 
//! for converting from the horzontal coordinate frame 
//! into the equatorial coordinate frame.

// import crate for matrix math
extern crate nalgebra as na;
extern crate astro;
extern crate chrono;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use chrono::{Datelike, Timelike};


#[derive(Debug, Clone)]
/// GPS struct
pub struct Gps {
    /// Latitude in degrees
    pub lat: f64,
    /// Longitude in degrees
    pub lon: f64,
    /// UTC time in Unix time (seconds from Jan 1 1970?)
    pub _utc: i64,
    /// greenwhich apparent sidereal time in degrees 
    /// (This is corrected for nutation and obliquity)
    pub gast: f64,
}

impl Gps {
    /// This function creates a new instance of the Gps struct
    /// initialized with selected initial values and then calculates
    /// the corresponding greenwhich apparent sidereal time (gast)
    pub fn new()-> Gps{
        let lat = 47.002282304922659;
        let lon = 82.210336016074365;
        let _utc = 946684800;
        let gast = 0.0;
        let mut gps = Gps{
            lat, 
            lon,
            _utc,
            gast,
        };
        gps.get_greenwhich_apparent_sidereal_time();
        gps
    }
    /// Calculates the greenwhich apparent mean sidereal time
    /// 
    /// # Detailed Explanation
    ///
    /// This functions updates the greenwhich apparent sidereal 
    /// time (GAST) while correcting for nutation and obliquity.
    /// The steps followed are:
    /// - get system time in Unix time (seconds since midnight Jan 1 1970)
    /// - calculate the current day in decimal form including the fraction 
    /// of a day that has passed
    /// - use decimal day in type astro::time::date to get a julian date (f64)
    /// - using julian date get mean sidereal time
    /// - using julian date get obliquity correction factor
    /// - using julian date get nutation correction factor
    /// - using the last 3 arguments get apparent sidereal time and convert to degrees
    /// 
    /// # Results
    /// - `self.gast: f64` - The greenwhich apparent sidereal time in degrees
    pub fn get_greenwhich_apparent_sidereal_time(&mut self){
        // this was tested against python and websites and thus is not in the matlab live script
        trace!("get_greenwhich_apparent_sidereal_time start");
        let dt = chrono::prelude::Utc::now();

        // break down the DateTime into a day object
        let _day = astro::time::DayOfMonth{
            day: dt.day() as u8,
            hr: dt.hour() as u8,
            min: dt.minute() as u8,
            sec: dt.second() as f64,
            time_zone: 0.0,
        };
        
        //use day object to calculate the decimal day
        let decimal_day = {
            _day.day as f64 
            + (_day.hr as f64)/24.0
            + _day.min as f64/(24.0*60.0)
            + _day.sec as f64/(24.0*3600.0)
        };

        // use decimal day to create a complete date object
        let _date = astro::time::Date {
            year: dt.year() as i16,
            month: dt.month() as u8,
            decimal_day,
            cal_type: astro::time::CalType::Gregorian,
        };
        
        //get the julian date
        let jd = astro::time::julian_day(&_date);

        //calculate the mean (uncorrected) apparent sidereal time
        let gmst = astro::time::mn_sidr(jd);

        // get the ecliptic oliquity correction factor
        let mn_obliq = astro::ecliptic::mn_oblq_laskar(jd);
        // get the nutation correction factor
        let nut = astro::nutation::nutation(jd);


        debug!("gmst {}, jd {}, mn_obliq {}, nut in lon {}, nut in obl {}", gmst.to_degrees(), jd, mn_obliq.to_degrees(), nut.0.to_degrees(), nut.1.to_degrees());

        let gast = astro::time::apprnt_sidr(gmst, nut.0, mn_obliq + nut.1);
        debug!("true obliq {}, apparent srt {}", (mn_obliq + nut.1).to_degrees(), gast.to_degrees());
        //convert to degrees and save in the struct
        self.gast = gast.to_degrees();
        trace!("get_greenwhich_apparent_sidereal_time end");
    }
}
