// import crate for matrix math
extern crate nalgebra as na;
extern crate astro;
extern crate chrono;

#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use chrono::{Datelike, Timelike};


#[derive(Debug, Clone)]
pub struct Gps {
    pub lat: f64,
    pub lon: f64,
    pub _utc: i64,
    pub gast: f64,
}

impl Gps {
    pub fn new()-> Gps{
        let lat = 47.002282304922659;
        let lon = 82.210336016074365;
        // let lat = 0.0;
        // let lon = 0.0;
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

    pub fn get_greenwhich_apparent_sidereal_time(&mut self){
        let dt = chrono::prelude::Utc::now();
        // let dt = chrono::DateTime::<chrono::Utc>::from_utc(
        //     chrono::NaiveDateTime::from_timestamp_opt(self.utc,0).unwrap(),
        //      chrono::Utc);
        let _day = astro::time::DayOfMonth{
            day: dt.day() as u8,
            hr: dt.hour() as u8,
            min: dt.minute() as u8,
            sec: dt.second() as f64,
            time_zone: 0.0,
        };
        
        let decimal_day = {
            _day.day as f64 
            + _day.hr as f64/24.0
            + _day.min as f64/(24.0*60.0)
            + _day.sec as f64/(24.0*3600.0)
        };

        let _date = astro::time::Date {
            year: dt.year() as i16,
            month: dt.month() as u8,
            decimal_day,
            cal_type: astro::time::CalType::Gregorian,
        };

        let jd = astro::time::julian_day(&_date);
        let gmst = astro::time::mn_sidr(jd);

        let mn_obliq = astro::ecliptic::mn_oblq_laskar(jd);


        let nut = astro::nutation::nutation(jd);


        debug!("gmst {}, jd {}, mn_obliq {}, nut in lon {}, nut in obl {}", gmst.to_degrees(), jd, mn_obliq.to_degrees(), nut.0.to_degrees(), nut.1.to_degrees());

        let gast = astro::time::apprnt_sidr(gmst, nut.0, mn_obliq + nut.1);
        debug!("true obliq {}, apparent srt {}", (mn_obliq + nut.1).to_degrees(), gast.to_degrees());
        self.gast = gast.to_degrees();
    }
}
