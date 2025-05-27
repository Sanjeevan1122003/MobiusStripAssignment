# Imported all the requried libraries
import numpy as np # Library for numerical operations
import matplotlib.pyplot as plt # Library for plotting
from scipy.integrate import simpson # Library for numerical integration


class MobiusStrip:
    # Class representing a Mobius strip in 3D space.
    def __init__(self, R=1.0, w=0.2, n=100):
        # Initialize the Mobius strip parameters.
        # Parameters:
        # R (float): Radius of the Möbius strip center circle.
        # w (float): Width of the strip.
        # n (int): Resolution for mesh generation.
        self.R = R
        self.w = w
        self.n = n
        self.u = np.linspace(0, 2 * np.pi, n)
        self.v = np.linspace(-w / 2, w / 2, n)
        self.U, self.V = np.meshgrid(self.u, self.v)
        self.X, self.Y, self.Z = self._compute_mesh()
        print("\nInitialized Mobius Strip")
        print(f"Radius = {self.R}")
        print(f"Width = {self.w}")
        print(f"Resolution = {self.n}")

    def _compute_mesh(self):
        # Compute the 3D mesh coordinates for the Mobius strip.
        # Returns:
        # Tuple of arrays (X, Y, Z) representing the surface.
        U, V = self.U, self.V
        X = (self.R + V * np.cos(U / 2)) * np.cos(U)
        Y = (self.R + V * np.cos(U / 2)) * np.sin(U)
        Z = V * np.sin(U / 2)
        return X, Y, Z
    
    # Function to approximate the surface area of the Mobius strip using numerical integration.
    # This Funtion
    # Returns:
    # float: Surface area.
    # Partial derivatives of the parametric surface
    def compute_surface_area(self):
        ru = np.empty((self.n, self.n, 3))
        rv = np.empty((self.n, self.n, 3))
        # ru: partial derivative with respect to u
        ru[:, :, 0] = -np.sin(self.U) * (self.R + self.V * np.cos(self.U / 2)) - \
                      0.5 * self.V * np.sin(self.U / 2) * np.cos(self.U)
        ru[:, :, 1] = np.cos(self.U) * (self.R + self.V * np.cos(self.U / 2)) - \
                      0.5 * self.V * np.sin(self.U / 2) * np.sin(self.U)
        ru[:, :, 2] = 0.5 * self.V * np.cos(self.U / 2)
        # rv: partial derivative with respect to v
        rv[:, :, 0] = np.cos(self.U / 2) * np.cos(self.U)
        rv[:, :, 1] = np.cos(self.U / 2) * np.sin(self.U)
        rv[:, :, 2] = np.sin(self.U / 2)
        # Surface area = ∫∫ ||ru × rv|| dudv
        cross = np.cross(ru, rv)
        norm = np.linalg.norm(cross, axis=2)
        area = simpson(simpson(norm, self.v), self.u)
        return area
    
    # Function to compute the length of the edge of the Mobius strip.
    # Returns:Edge length.
    def compute_edge_length(self):
        edge_u = self.u
        v_edge = self.w / 2  # Top edge of the strip
        # Parametric equations of the edge
        x = (self.R + v_edge * np.cos(edge_u / 2)) * np.cos(edge_u)
        y = (self.R + v_edge * np.cos(edge_u / 2)) * np.sin(edge_u)
        z = v_edge * np.sin(edge_u / 2)
        # Compute differential arc lengths
        dx = np.gradient(x)
        dy = np.gradient(y)
        dz = np.gradient(z)
        integrand = np.sqrt(dx**2 + dy**2 + dz**2)
        length = simpson(integrand, edge_u)
        return length
    
    # Function to Plot the Möbius strip in 3D
    def plot(self, colormap='plasma'):
        # Parameters:
        # colormap (str): Matplotlib colormap for shading the surface.
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        # Color the surface based on Z-values
        norm = plt.Normalize(self.Z.min(), self.Z.max())
        colors = plt.get_cmap(colormap)(norm(self.Z))
        ax.plot_surface(self.X, self.Y, self.Z, facecolors=colors, edgecolor='none', alpha=1.0)
        # Axes labels
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        # Reference axes
        ax.quiver(0, 0, 0, 1, 0, 0, color='red', arrow_length_ratio=0.1)
        ax.quiver(0, 0, 0, 0, 1, 0, color='green', arrow_length_ratio=0.1)
        ax.quiver(0, 0, 0, 0, 0, 1, color='blue', arrow_length_ratio=0.1)
        # Annotation box with strip info
        info_text = (
            f"Radius = {self.R:.2f}\n"
            f"Width = {self.w:.2f}\n"
            f"Resolution = {self.n}\n"
            f"Surface Area = {self.compute_surface_area():.4f}\n"
            f"Edge Length = {self.compute_edge_length():.4f}"
        )
        fig.text(0.8, 0.95, info_text, ha='left', va='top', fontsize=10,
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

        plt.tight_layout()
        plt.show()


def main():
    # Main execution function for interactive input and visualization.
    print("Please enter only numreic values for radius, width, and resolution.")
    try:
        R_input = input("\nEnter the Radius: ")
        w_input = input("Enter the Width: ")
        n_input = input("Enter the Resolution: ")

        R = float(R_input)
        w = float(w_input)
        n = int(n_input)
    except ValueError:
        print("Invalid input. Please enter numeric values.")
        return
    # Create a Mobius strip object with the given input parameters
    mobius = MobiusStrip(R, w, n)
    # Outputs of the class methods
    print("\nOutputs:")
    print(f"Surface Area = {mobius.compute_surface_area():.4f}")
    print(f"Edge Length = {mobius.compute_edge_length():.4f}")
    mobius.plot(colormap='Blues_r')

# Running the main function
if __name__ == "__main__":
    main()
