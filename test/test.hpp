namespace gyper
{

class Logger
{
//public:
//  static Options * instance(); // Gets the only instance (singleton pattern)
//  static const Options * const_instance(); // const version of instance()

private:
  Logger(); // Prevent construction of new instances
  Logger(Logger const &); // Prevent copy-construction
  Logger & operator=(Logger const &); // Prevent copy-assignment

  static Logger * _instance;
};

} // namespace gyper
