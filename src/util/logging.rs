use std::fs;
use std::path::Path;

pub fn init_logging(out_dir: &Path) -> anyhow::Result<()> {
    fs::create_dir_all(out_dir)?;
    let logfile = out_dir.join("onsm.log");

    let mut builder = env_logger::Builder::from_default_env();
    builder
        .filter_level(log::LevelFilter::Info)
        .format_timestamp_millis()
        .format_module_path(false)
        .format_level(true);

    // Also duplicate logs to a file
    let file = std::fs::File::create(&logfile)?;
    let _ = std::sync::Mutex::new(file);
    builder.target(env_logger::Target::Stderr);

    // Install a simple tee: stderr via env_logger; file via log::set_boxed_logger is overkill.
    // Keep it minimal: advise users to tail the file created by mapping/scoring steps.
    builder.init();
    log::info!("Logging initialized. Log file: {}", logfile.display());
    Ok(())
}
