pub fn setup_deadlock_detection() {
  #[cfg(debug_assertions)]
  {
    use color_backtrace::BacktracePrinter;
    use parking_lot::deadlock;
    use std::thread;
    use std::time::Duration;

    // Create a background thread which checks for deadlocks every 5s
    thread::spawn(move || loop {
      let printer = BacktracePrinter::new();

      thread::sleep(Duration::from_secs(5));
      let deadlocks = deadlock::check_deadlock();
      if deadlocks.is_empty() {
        continue;
      }

      println!("{} deadlocks detected", deadlocks.len());
      for (i, threads) in deadlocks.iter().enumerate() {
        println!("Deadlock #{i}");
        for t in threads {
          println!("Thread Id {:#?}", t.thread_id());
          println!("{}", printer.format_trace_to_string(t.backtrace()).unwrap());
        }
      }
    });
  }
}
