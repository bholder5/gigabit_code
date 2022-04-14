use chrono::Local;
use env_logger::Builder;
pub use log::LevelFilter;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
pub use log4rs::append::file::FileAppender;
pub use log4rs::config::{Appender, Config, Root};
pub use log4rs::encode::pattern::PatternEncoder;
use std::io::Write;

// Initiates a log file and format
// initiate log/output.log as the output file for our log.
pub fn init_log() -> () {
    let logfile = FileAppender::builder()
        .encoder(Box::new(PatternEncoder::new("{l} - {m}\n")))
        .build("log/output.log");

    let file_avail = match &logfile {
        Ok(_file) => 1,
        Err(_error) => 0,
    };

    //flags for error checking the config function and init functions work for logger
    //if not, env_logger will be activated
    let mut config_avail = 0;
    let mut init_avail = 0;
    let mut fail = false; //flag used for initializing environmental logger if file initialization fails

    if file_avail == 1 {
        let logfile = match logfile {
            Ok(file) => file,
            Err(error) => panic!(
                "This will never happen but need to do it to make the compiler happy. {}",
                error
            ),
        };

        let config = Config::builder()
            .appender(Appender::builder().build("logfile", Box::new(logfile)))
            .build(
                Root::builder()
                    .appender("logfile")
                    .build(LevelFilter::Off),
            );

        // unpack config and set flag based on success (1) or not (0)
        config_avail = match &config {
            Ok(_file) => 1,
            Err(_error) => 0,
        };

        // unpack config
        if config_avail == 1 {
            let config = match config {
                Ok(config) => config,
                //never happens because config is confirmed prior to this
                Err(error) => panic!(
                    "This will never happen but need to do it to make the compiler happy. {}",
                    error
                ),
            };

            let log_init = log4rs::init_config(config);

            init_avail = match &log_init {
                Ok(_handle) => 1,
                Err(_error) => 0,
            };

            if init_avail == 1 {
                let _log_init = match log_init {
                    Ok(handle) => handle,
                    Err(error) => panic!(
                        "This will never happen but need to do it to make the compiler happy. {}",
                        error
                    ),
                };
            } else {
                error!("Log file initialization failed");
                fail = true;
            }
        } else {
            error!("Log file Config failed");
            fail = true;
        }
    } else {
        fail = true;
    }

    if fail == true {
        error!("Could not open the log file for appending");

        //initialize an evironment logger
        Builder::new()
            .format(|buf, record| {
                writeln!(
                    buf,
                    "{} [{}] - {}",
                    Local::now().format("%Y-%m-%dT%H:%M:%S"),
                    record.level(),
                    record.args()
                )
            })
            .filter(None, LevelFilter::Off)
            .init();
    }

    info!("\n***************************\n***************************\n NEW RUN STARTS HERE \n*************************** \n***************************");
    info!("File logging success flags: \n file available: {}, \n config available: {}, \n init available: {}", file_avail, config_avail, init_avail);
}
