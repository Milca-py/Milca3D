import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QOpenGLWidget
from PyQt5.QtCore import Qt
from OpenGL.GL import *
from OpenGL.GLU import *

class StructuralViewer(QOpenGLWidget):
    def initializeGL(self):
        glEnable(GL_DEPTH_TEST)
        glClearColor(0.9, 0.9, 0.9, 1)

    def resizeGL(self, w, h):
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, w / h if h else 1, 0.1, 1000.0)
        glMatrixMode(GL_MODELVIEW)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        gluLookAt(5, 5, 5, 0, 0, 0, 0, 1, 0)

        # Modelo estructural básico (ejes de referencia)
        glBegin(GL_LINES)
        # Eje X (rojo)
        glColor3f(1, 0, 0)
        glVertex3f(0, 0, 0)
        glVertex3f(3, 0, 0)
        # Eje Y (verde)
        glColor3f(0, 1, 0)
        glVertex3f(0, 0, 0)
        glVertex3f(0, 3, 0)
        # Eje Z (azul)
        glColor3f(0, 0, 1)
        glVertex3f(0, 0, 0)
        glVertex3f(0, 0, 3)
        glEnd()

        # View Cube
        self.draw_view_cube()

    def draw_view_cube(self):
        glPushMatrix()
        glTranslatef(2.8, 2.8, -5.5)  # Posición relativa
        glScalef(0.3, 0.3, 0.3)
        glBegin(GL_QUADS)

        faces = [
            ([1, 0, 0], [1, 1, 1], [1, -1, 1], [1, -1, -1], [1, 1, -1]),  # X+
            ([0, 1, 0], [-1, 1, 1], [1, 1, 1], [1, 1, -1], [-1, 1, -1]),  # Y+
            ([0, 0, 1], [-1, 1, 1], [1, 1, 1], [1, -1, 1], [-1, -1, 1]),  # Z+
        ]

        for color, v1, v2, v3, v4 in faces:
            glColor3fv(color)
            glVertex3fv(v1)
            glVertex3fv(v2)
            glVertex3fv(v3)
            glVertex3fv(v4)

        glEnd()
        glPopMatrix()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Visualizador de Modelo Estructural")
        self.setGeometry(100, 100, 800, 600)
        self.viewer = StructuralViewer()
        self.setCentralWidget(self.viewer)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())
