#include "includes.h"

extern GLubyte *mybuffer;

void display()
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glDrawPixels(width, height, GL_RGB, GL_UNSIGNED_BYTE, mybuffer);

    glutSwapBuffers();
    glutPostRedisplay();
}

int main(int argc, char **argv)
{
    int WindowID;

    mybuffer = (GLubyte *)malloc(width * height * 3 * sizeof(GLubyte));

    for (int i = 0; i < width * height * 3; i++)
    {
        mybuffer[i]     = 0;
        mybuffer[i + 1] = 0;
        mybuffer[i + 2] = 0;
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(width, height);
    WindowID = glutCreateWindow("Ray Tracer");

    MakeGUI();

    glutDisplayFunc(display);
    glutMainLoop();

    free(mybuffer);
    return 0;
}