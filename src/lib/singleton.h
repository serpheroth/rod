#ifndef SINGLETON_H
#define SINGLETON_H

template <class T>
class Singleton
{
public:
       static inline T* Instance() {
         static T instance;
         return &instance;
       }
private:
       Singleton(void){}
       ~Singleton(void){}
       Singleton(const Singleton&){}
       Singleton & operator= (const Singleton &){}
};

//Class that will implement the singleton mode,
//must use the macro in it's delare file
#define DECLARE_SINGLETON_CLASS(type) \
       friend class Singleton<type>; \

#endif // SINGLETON_H
