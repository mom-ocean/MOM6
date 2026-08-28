.. autodoxymodule:: {{ fullname }}
   :members:
   {% if methods %}
   :methods:
   {% endif %}
   {% if types %}
   :types:
   {% endif %}

   {% if types %}
   ----------
   Data Types
   ----------

   .. autodoxysummary::
      :kind: type

   {% for item in types %}
      ~{{ item }}
   {% endfor %}
   {% endif %}

   {% if methods %}
   ---------------------
   Functions/Subroutines
   ---------------------

   .. autodoxysummary::
      :kind: func

   {% for item in methods %}
      ~{{ fullname }}::{{ item }}
   {% endfor %}
   {% endif %}
