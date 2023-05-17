# Servidor

---

## Puesta en marcha del servidor

Si es usted el que hostea la pagina y debe levantar el servidor Genie, debe realizar los siguientes pasos:

1.  Diríjase a la carpeta bin la cual estará dentro de la carpeta de la pagina web, sabrá que a llegado a ella porque encontrara los archivos repl, server ... etc.
2.  Ejecute el siguiente comando, que es para iniciar la aplicación y salir de la consola sin cerrar la aplicación:

        screen -S genie

3.  Luego establezca la GENIE_ENV variable de entorno en prod:

        export GENIE_ENV=prod 
        
4.  Ejecute el archvo server:

        ./server
        
5.  Salga de screen con:

        CTRL+d 

screen -X -S genie quit, esto mata el proceso, CTRL+a CTRL+d, CTRL+a d
